# dod.py
# Author: Jagger Hershey
#
# All analysis specific to the DoD/VA chemicals-of-concern list lives in this one script, kept
# in clearly labeled steps/sections so it's obvious what each part does and depends on:
#
#   STEP 1 (done): match the DoD list to ChEMBL compounds -> dod.spl, dod.chembl_match,
#                  dod.chembl_activities (bioactivity + mechanism-of-action data)
#   STEP 2 (done): commonalities across the DoD list in CTD -> shared/specific genes,
#                  shared diseases (independent of whether a chemical matched ChEMBL - CTD name
#                  coverage and ChEMBL compound coverage are separate concerns)
#   STEP 3 (done): BindingDB binding affinity per ChEMBL-matched DoD chemical (decoupled from
#                  CTD hub-gene status - restricting to hub genes returned zero matches, a real
#                  finding, not a bug: see binding_affinity_per_chemical's docstring)
#   STEP 4 (done): per-chemical gene-disease NetworkX graphs, saved to data/dod/networks/
#   STEP 5 (done): FDA-approved counteracting-drug candidates from ChEMBL, for every gene a
#                  chemical affects (not just the top-20 subset drawn in step 4's networks)
#   STEP 6 (done): one aggregated HTML report - a summary section (commonalities across the
#                  list) plus one section per chemical (ChEMBL/BindingDB info, gene evidence,
#                  the network image, and top counteracting drugs), saved to data/dod/dod_report.html
#
# This is a sibling of analyze.py/report.py (the ATSDR-based CTD analysis), not a replacement -
# it writes to its own `dod` schema and never touches atsdr.*.
#
# ChEMBL itself is never copied in - it's queried directly from its own SQLite file via the
# `chembl` schema attached read-only by scripts/common/db.py.
import re
import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'common'))
from db import connect, PROJECT_ROOT

DOD_CSV = PROJECT_ROOT / 'data' / 'dod' / 'DoD-Chemicals-of-Concern.csv'
NETWORKS_OUT_DIR = PROJECT_ROOT / 'data' / 'dod' / 'networks'
PUBCHEM_CACHE_PATH = PROJECT_ROOT / 'data' / 'chembl' / 'raw' / 'pubchem_casrn_inchikey_cache.json'

PUBCHEM_CID_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/{casrn}/cids/JSON"
PUBCHEM_INCHIKEY_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON"
PUBCHEM_REQUEST_DELAY_S = 0.25  # stay under PubChem's ~5 requests/second guidance

# =========================================================================================
# STEP 1: Match DoD chemicals to ChEMBL compounds (bioactivity + mechanism-of-action data)
# =========================================================================================

# Bioactivity + mechanism-of-action data for every matched molecule in one shot, joined across:
# compound structure -> activities -> assay -> target -> target's protein component (UniProt)
# -> curated mechanism-of-action (drug_mechanism), where one exists for that compound/target pair.
ACTIVITY_JOIN_SQL = """
CREATE OR REPLACE TABLE dod.chembl_activities AS
SELECT
    mm.dod_substance_name,
    mm.dod_casrn,
    mm.match_method,
    md.chembl_id            AS compound_chembl_id,
    md.pref_name             AS compound_pref_name,
    cs.canonical_smiles,
    cs.standard_inchi_key,
    act.activity_id,
    act.standard_type,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.pchembl_value,
    act.action_type          AS activity_action_type,
    a.assay_type,
    a.description             AS assay_description,
    a.confidence_score,
    td.chembl_id               AS target_chembl_id,
    td.pref_name                AS target_pref_name,
    td.target_type,
    td.organism                  AS target_organism,
    cseq.accession                 AS target_uniprot_accession,
    cseq.db_source                  AS target_seq_db_source,
    dm.mechanism_of_action,
    dm.action_type            AS mechanism_action_type,
    dm.binding_site_comment,
    dm.mechanism_comment,
    dm.direct_interaction,
    dm.molecular_mechanism,
    doc.pubmed_id,
    doc.doi
FROM _matched_molregnos mm
JOIN chembl.molecule_dictionary md ON md.molregno = mm.molregno
JOIN chembl.compound_structures cs ON cs.molregno = md.molregno
JOIN chembl.activities act ON act.molregno = md.molregno
JOIN chembl.assays a ON a.assay_id = act.assay_id
LEFT JOIN chembl.target_dictionary td ON td.tid = a.tid
LEFT JOIN chembl.target_components tc ON tc.tid = td.tid
LEFT JOIN chembl.component_sequences cseq ON cseq.component_id = tc.component_id
LEFT JOIN chembl.drug_mechanism dm ON dm.molregno = md.molregno AND dm.tid = td.tid
LEFT JOIN chembl.docs doc ON doc.doc_id = a.doc_id
"""


def _normalize_name(name: str) -> str:
    return re.sub(r'[^A-Z0-9]', '', name.upper())


def load_dod_list(con):
    con.execute(f"""
        CREATE OR REPLACE TABLE dod.spl AS
        SELECT * FROM read_csv('{DOD_CSV.as_posix()}', header=true)
    """)
    return con.execute('SELECT "Substance Name", "CASRN" FROM dod.spl').fetchall()


def build_cas_and_name_indexes(con, dod_list):
    """Single pass over molecule_synonyms/molecule_dictionary (via the attached chembl schema),
    matched against the (small) DoD list, instead of one query per substance against
    multi-million-row tables.

    Keyed by substance NAME rather than CASRN - unlike ATSDR, several DoD entries (Fluorocarbons,
    Particulate Matter, ...) are categories with no single CASRN, so CASRN can't be trusted as a
    unique key here (multiple blank-CASRN substances would otherwise collide on the same key)."""
    casdigits_to_name = {re.sub(r'\D', '', casrn): name for name, casrn in dod_list if casrn}
    norm_name_to_name = {_normalize_name(name): name for name, _ in dod_list if name}

    cas_matches = {}
    name_matches = {}

    for molregno, synonym in con.execute("SELECT molregno, synonyms FROM chembl.molecule_synonyms WHERE synonyms IS NOT NULL").fetchall():
        digits = re.sub(r'\D', '', synonym)
        if digits and digits in casdigits_to_name:
            cas_matches.setdefault(casdigits_to_name[digits], molregno)
        norm = _normalize_name(synonym)
        if norm and norm in norm_name_to_name:
            name_matches.setdefault(norm_name_to_name[norm], molregno)

    for molregno, pref_name in con.execute("SELECT molregno, pref_name FROM chembl.molecule_dictionary WHERE pref_name IS NOT NULL").fetchall():
        norm = _normalize_name(pref_name)
        if norm and norm in norm_name_to_name:
            name_matches.setdefault(norm_name_to_name[norm], molregno)

    return cas_matches, name_matches


def _pubchem_inchikey_for_casrn(casrn: str, cache: dict):
    if casrn in cache:
        return cache[casrn]
    inchikey = None
    try:
        url = PUBCHEM_CID_URL.format(casrn=urllib.parse.quote(casrn))
        with urllib.request.urlopen(url, timeout=15) as resp:
            cid = json.load(resp)['IdentifierList']['CID'][0]
        time.sleep(PUBCHEM_REQUEST_DELAY_S)
        url = PUBCHEM_INCHIKEY_URL.format(cid=cid)
        with urllib.request.urlopen(url, timeout=15) as resp:
            inchikey = json.load(resp)['PropertyTable']['Properties'][0]['InChIKey']
        time.sleep(PUBCHEM_REQUEST_DELAY_S)
    except Exception as exc:
        print(f"  PubChem lookup failed for CASRN {casrn}: {exc}")
        inchikey = None
    cache[casrn] = inchikey
    return inchikey


def match_remaining_via_pubchem(con, dod_list, already_matched: set):
    """Returns a dict keyed by substance NAME (not CASRN - see build_cas_and_name_indexes for
    why). Substances with a blank CASRN are skipped here since there's nothing to look up."""
    remaining = [(name, casrn) for name, casrn in dod_list if casrn and name not in already_matched]
    if not remaining:
        return {}

    cache = json.load(open(PUBCHEM_CACHE_PATH)) if PUBCHEM_CACHE_PATH.exists() else {}

    print(f"Resolving {len(remaining)} unmatched substances via PubChem CASRN -> InChIKey ...")
    inchikey_to_name = {}
    for name, casrn in remaining:
        inchikey = _pubchem_inchikey_for_casrn(casrn, cache)
        if inchikey:
            inchikey_to_name[inchikey] = name
    with open(PUBCHEM_CACHE_PATH, 'w') as f:
        json.dump(cache, f, indent=2)

    if not inchikey_to_name:
        return {}

    matches = {}
    for molregno, inchikey in con.execute("SELECT molregno, standard_inchi_key FROM chembl.compound_structures WHERE standard_inchi_key IS NOT NULL").fetchall():
        if inchikey in inchikey_to_name:
            matches.setdefault(inchikey_to_name[inchikey], molregno)
    return matches


def match_dod_chemicals(con):
    """Step 1: load the DoD list and match each substance to a ChEMBL compound, writing
    dod.spl, dod.chembl_match (every substance, matched or not) and dod.chembl_activities
    (bioactivity + mechanism-of-action data for whatever matched).

    matches is keyed by substance NAME throughout (not CASRN) - several DoD entries share a
    blank CASRN (they're categories like Fluorocarbons/Particulate Matter, not single defined
    compounds), and name is the only field guaranteed unique across this list."""
    dod_list = [(name.strip(), casrn.strip() if casrn else '') for name, casrn in load_dod_list(con)]
    casrn_by_name = {name: casrn for name, casrn in dod_list}

    cas_matches, name_matches = build_cas_and_name_indexes(con, dod_list)
    matches = {name: (molregno, 'cas_synonym') for name, molregno in cas_matches.items()}
    for name, molregno in name_matches.items():
        matches.setdefault(name, (molregno, 'name_match'))

    pubchem_matches = match_remaining_via_pubchem(con, dod_list, set(matches))
    for name, molregno in pubchem_matches.items():
        matches.setdefault(name, (molregno, 'pubchem_inchikey'))

    # dod.chembl_match: every DoD substance, matched or not
    con.execute("""
        CREATE OR REPLACE TABLE dod.chembl_match (
            substance_name VARCHAR, casrn VARCHAR, chembl_molregno BIGINT,
            match_method VARCHAR, matched BOOLEAN
        )
    """)
    match_rows = [
        (name, casrn, matches[name][0] if name in matches else None,
         matches[name][1] if name in matches else None, name in matches)
        for name, casrn in dod_list
    ]
    con.executemany(
        "INSERT INTO dod.chembl_match VALUES (?, ?, ?, ?, ?)", match_rows
    )

    # One SQL join across the attached chembl schema for every matched compound at once,
    # instead of looping per-molregno in Python.
    con.execute("CREATE OR REPLACE TEMP TABLE _matched_molregnos (molregno BIGINT, dod_substance_name VARCHAR, dod_casrn VARCHAR, match_method VARCHAR)")
    con.executemany(
        "INSERT INTO _matched_molregnos VALUES (?, ?, ?, ?)",
        [(molregno, name, casrn_by_name[name], method) for name, (molregno, method) in matches.items()],
    )
    con.execute(ACTIVITY_JOIN_SQL)

    total_rows = con.execute("SELECT count(*) FROM dod.chembl_activities").fetchone()[0]
    substances_with_activity = con.execute("SELECT count(DISTINCT dod_substance_name) FROM dod.chembl_activities").fetchone()[0]
    method_counts = dict(con.execute("SELECT match_method, count(*) FROM dod.chembl_match WHERE matched GROUP BY match_method").fetchall())

    print(f"Done: dod.chembl_activities, dod.chembl_match")
    print(f"Matched {len(matches)}/{len(dod_list)} DoD substances to a ChEMBL compound")
    print(f"  of which {substances_with_activity} had activity/bioactivity data")
    print(f"Match methods: {method_counts}")
    print(f"Total activity rows written: {total_rows:,}")


# =========================================================================================
# STEP 2: Commonalities across the DoD list in CTD (shared/specific genes, shared diseases)
# =========================================================================================
# Mirrors analyze.py's atsdr_top_chemical_patterns() gene/disease-sharing logic, but scoped to
# the DoD list instead of ATSDR's ranked top-N (the DoD list has no Rank column - it's already
# a small, curated set, so there's nothing to subset further). Uses the whole 35-chemical list
# from dod.spl, not just the ChEMBL-matched subset from Step 1 - CTD name coverage and ChEMBL
# compound coverage are independent concerns, and CTD's chem_gene_ixns/chemicals_diseases
# already have exact-name matches for all 35
#
# Requires dod.spl to exist - i.e. match_dod_chemicals(con) (or at least load_dod_list(con))
# must have already run this session.

def _match_dod_to_ctd(con, ctd_table: str):
    """Match every DoD substance to an exact (case-insensitive) ChemicalName in the given CTD
    table (chem_gene_ixns or chemicals_diseases), reporting any that don't resolve as a coverage
    gap instead of silently dropping them - a substance could exist in CTD under a different
    exact name (e.g. a specific isomer) that this won't catch, and that's worth knowing about
    rather than hiding. Populates temp table _dod_ctd_matched(Substance, ChemicalName)."""
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _dod_ctd_matched AS
        SELECT DISTINCT d."Substance Name" AS Substance, t.ChemicalName
        FROM dod.spl d JOIN ctd.{ctd_table} t ON lower(t.ChemicalName) = lower(d."Substance Name")
    """)
    matched = con.execute("SELECT DISTINCT Substance FROM _dod_ctd_matched").df()
    all_dod = con.execute('SELECT "Substance Name" AS Substance FROM dod.spl').df()
    unmatched = all_dod[~all_dod['Substance'].isin(matched['Substance'])]
    return matched, unmatched


def common_genes(con, human_only=True, min_shared=10, display_n=20):
    """
    Within the DoD chemical list, find:
      - genes hit by >= min_shared of the DoD chemicals ("hub" genes within this list), broken
        down by which interaction type (Suffix) each contributing chemical hits that gene through
      - genes hit by exactly one DoD chemical ("specific" to that one chemical)

    min_shared defaults to 10 rather than 2 - at 2, CTD's heavy literature coverage of
    well-studied chemicals (Lead, Cadmium, ...) means almost any two substances share
    thousands of genes incidentally, so "hub" ends up meaning almost nothing (it split
    ~15,800 genes roughly in half). 10 is still just a fixed cutoff, not derived from the
    actual distribution - revisit if it stops feeling selective enough.
    """
    label = "Human-Only" if human_only else "All Organisms"
    organism_filter = "AND Organism ILIKE '%Homo sapiens%'" if human_only else ""

    matched, unmatched = _match_dod_to_ctd(con, 'chem_gene_ixns')
    print(f"===== DoD List -> CTD chem_gene_ixns Chemical Match ({label}) =====")
    print(f"Matched {len(matched)}/{len(matched) + len(unmatched)} to a CTD chemical name")
    if len(unmatched):
        print("Unmatched (no exact ChemicalName match in ctd.chem_gene_ixns - coverage gap, not analyzed further):")
        print(unmatched)

    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _dod_chem_gene_split AS
        SELECT dm.Substance, cg.GeneSymbol,
               lower(trim(split_part(UNNEST(string_split(cg.InteractionActions, '|')), '^', 2))) AS Suffix
        FROM ctd.chem_gene_ixns cg
        JOIN _dod_ctd_matched dm ON dm.ChemicalName = cg.ChemicalName
        WHERE cg.InteractionActions IS NOT NULL {organism_filter}
    """)
    gene_share = con.execute("""
        SELECT Substance, GeneSymbol, Suffix, count(*) AS n
        FROM _dod_chem_gene_split GROUP BY Substance, GeneSymbol, Suffix
    """).df()

    n_chemicals_per_gene = gene_share.groupby('GeneSymbol')['Substance'].nunique().rename('n_chemicals')

    hub_genes = n_chemicals_per_gene[n_chemicals_per_gene >= min_shared].sort_values(ascending=False)
    print(f"\n===== Genes Shared Across >= {min_shared} DoD Chemicals ({label}) =====")
    print(f"Genes hit by >= {min_shared} of the {len(matched)} matched DoD chemicals: {len(hub_genes)}")
    for gene, n_chem in hub_genes.head(display_n).items():
        print(f"\n{gene} (hit by {n_chem} chemicals):")
        gene_rows = gene_share[gene_share['GeneSymbol'] == gene]
        for substance, grp in gene_rows.groupby('Substance'):
            types = ', '.join(sorted(grp['Suffix'].unique()))
            print(f"    {substance}: {types}")

    specific_genes = n_chemicals_per_gene[n_chemicals_per_gene == 1]
    specific_detail = gene_share[gene_share['GeneSymbol'].isin(specific_genes.index)][['Substance', 'GeneSymbol']].drop_duplicates()
    per_chem_specific_counts = specific_detail.groupby('Substance').size().sort_values(ascending=False)

    print(f"\n===== Genes Specific to a Single DoD Chemical ({label}) =====")
    print(f"Genes hit by exactly 1 of the matched DoD chemicals: {len(specific_genes)}")
    print("Specific-gene count per chemical (how many genes are touched by that chemical alone):")
    print(per_chem_specific_counts.head(display_n))

    return {"matched": matched, "unmatched": unmatched, "hub_genes": hub_genes, "specific_genes": specific_genes}


def common_diseases(con, min_shared=2, display_n=20):
    """Diseases linked (DirectEvidence only) to >= min_shared of the DoD chemicals in
    ctd.chemicals_diseases."""
    matched, unmatched = _match_dod_to_ctd(con, 'chemicals_diseases')
    print(f"===== DoD List -> CTD chemicals_diseases Chemical Match =====")
    print(f"Matched {len(matched)}/{len(matched) + len(unmatched)} to a CTD chemical name")
    if len(unmatched):
        print("Unmatched (no exact ChemicalName match in ctd.chemicals_diseases - coverage gap, not analyzed further):")
        print(unmatched)

    # "Death" is excluded - CTD records it as a disease term, but it's a non-specific outcome
    # rather than a clinically meaningful diagnosis, and shows up as a top "shared disease" for
    # almost any toxic chemical without saying anything useful about mechanism.
    disease_pairs = con.execute("""
        SELECT DISTINCT dm.Substance, cd.DiseaseName
        FROM ctd.chemicals_diseases cd
        JOIN _dod_ctd_matched dm ON dm.ChemicalName = cd.ChemicalName
        WHERE cd.DirectEvidence IS NOT NULL AND cd.DiseaseName != 'Death'
    """).df()
    n_chemicals_per_disease = disease_pairs.groupby('DiseaseName')['Substance'].nunique().rename('n_chemicals')
    shared_diseases = n_chemicals_per_disease[n_chemicals_per_disease >= min_shared].sort_values(ascending=False)

    print(f"\n===== Diseases Shared Across >= {min_shared} DoD Chemicals =====")
    print(f"Diseases linked (direct evidence) to >= {min_shared} of the {len(matched)} matched DoD chemicals: {len(shared_diseases)}")
    for disease, n_chem in shared_diseases.head(display_n).items():
        subs = sorted(disease_pairs.loc[disease_pairs['DiseaseName'] == disease, 'Substance'].unique())
        print(f"  {disease} ({n_chem}): {', '.join(subs)}")

    return {"matched": matched, "unmatched": unmatched, "shared_diseases": shared_diseases, "disease_pairs": disease_pairs}


# =========================================================================================
# STEP 3: BindingDB binding affinity for the DoD chemicals
# =========================================================================================
# Originally scoped to only the (chemical, hub gene) pairs from step 2, restricted to genes
# shared by >= min_shared DoD chemicals in CTD. That came back with ZERO matches - not a bug
# (verified both halves of the join independently): CTD's hub genes are broad toxicology-
# literature associations (TNF, CXCL8, IL6, CASP3, TP53, ...), while BindingDB only records
# direct, quantitative, competitive binding assays, which for these chemicals land on classic
# pharmacology targets instead (GPCRs, carbonic anhydrase, transporters, acetylcholinesterase
# for the nerve agents). The two are genuinely different slices of biology for toxicants (as
# opposed to drugs), and that non-overlap is itself worth knowing - see the overlap check below.
#
# So this reports ALL BindingDB binding data per DoD chemical, decoupled from CTD hub-gene
# status, since that's real per-chemical mechanism data worth keeping regardless. Ligand
# identity is matched via InChIKey rather than ChEMBL ID - BindingDB's own "ChEMBL ID of
# Ligand" field is only ~35% filled, so it undercounts real coverage (confirmed: ChEMBL ID
# matching found data for 2/21 chemicals, InChIKey matching found 5/21).
#
# TWO independent InChIKey sources are unioned, because they don't agree with each other:
#   - chembl_inchikey: dod.chembl_match's molregno -> chembl.compound_structures.standard_inchi_key
#   - ctd_inchikey: ctd.chemicals.InChIKey (CTD's own vocabulary file) matched directly against
#     BindingDB, bypassing ChEMBL entirely
# Confirmed empirically that these disagree even for compounds both have matched: CTD and ChEMBL
# sometimes compute different InChIKeys for "the same" chemical (different ionization/salt form
# chosen as canonical - seen concretely for Lead, Cadmium, Chromium, Copper, Silver). Using only
# one source silently drops real coverage the other would have found (e.g. Tabun and Toluene
# 2,4-Diisocyanate only turn up via the ctd_inchikey route). SMILES was also tested as a
# candidate match key and performed WORSE than either InChIKey route (same top 2 hits, lost
# Sarin/Soman) - BindingDB and ChEMBL evidently canonicalize SMILES with different software, so
# it's not used here.
#
# Every DoD substance is reported with which route(s) found it (or neither) - the "neither"
# list is for manual follow-up (e.g. looking the chemical up on BindingDB's site directly),
# since a blank result here could still be a matching/identifier gap rather than true absence.
#
# Requires match_dod_chemicals(con) to have already run this session (dod.chembl_match and
# dod.spl). The CTD hub-gene overlap check additionally requires common_genes(con) to have run
# (_dod_chem_gene_split).

def binding_affinity_per_chemical(con, display_n=20):
    """For each DoD chemical, pull ALL BindingDB binding affinity records (Ki/IC50/Kd/EC50)
    against any target, via the union of two independent InChIKey routes. Writes
    dod.binding_affinity (tagged with match_route), prints per-chemical route coverage (for
    manual follow-up on anything found via neither route) and a sample of records, then checks
    how many (if any) land on a gene CTD also reports an interaction with for that chemical."""
    con.execute("""
        CREATE OR REPLACE TABLE dod.binding_affinity AS
        WITH chembl_route AS (
            SELECT
                cm.substance_name AS dod_substance_name,
                md.chembl_id AS compound_chembl_id,
                bdb."Target Name" AS target_name,
                bdb."UniProt (SwissProt) Primary ID of Target Chain 1" AS uniprot_swissprot,
                bdb."UniProt (TrEMBL) Primary ID of Target Chain 1" AS uniprot_trembl,
                bdb."Ki (nM)" AS ki_nm,
                bdb."IC50 (nM)" AS ic50_nm,
                bdb."Kd (nM)" AS kd_nm,
                bdb."EC50 (nM)" AS ec50_nm,
                bdb."Curation/DataSource" AS source,
                bdb."PMID" AS pmid,
                'chembl_inchikey' AS match_route
            FROM dod.chembl_match cm
            JOIN chembl.molecule_dictionary md ON md.molregno = cm.chembl_molregno
            JOIN chembl.compound_structures cs ON cs.molregno = md.molregno
            JOIN bindingdb.activities bdb ON bdb."Ligand InChI Key" = cs.standard_inchi_key
            WHERE cm.matched
        ),
        ctd_route AS (
            SELECT
                d."Substance Name" AS dod_substance_name,
                NULL AS compound_chembl_id,
                bdb."Target Name" AS target_name,
                bdb."UniProt (SwissProt) Primary ID of Target Chain 1" AS uniprot_swissprot,
                bdb."UniProt (TrEMBL) Primary ID of Target Chain 1" AS uniprot_trembl,
                bdb."Ki (nM)" AS ki_nm,
                bdb."IC50 (nM)" AS ic50_nm,
                bdb."Kd (nM)" AS kd_nm,
                bdb."EC50 (nM)" AS ec50_nm,
                bdb."Curation/DataSource" AS source,
                bdb."PMID" AS pmid,
                'ctd_inchikey' AS match_route
            FROM dod.spl d
            JOIN ctd.chemicals cc ON lower(cc.ChemicalName) = lower(d."Substance Name")
            JOIN bindingdb.activities bdb ON bdb."Ligand InChI Key" = cc.InChIKey
            WHERE cc.InChIKey IS NOT NULL
        )
        SELECT * FROM chembl_route
        UNION
        SELECT * FROM ctd_route
    """)

    total_rows = con.execute("SELECT count(*) FROM dod.binding_affinity").fetchone()[0]

    # Per-substance route coverage: found via chembl route, ctd route, both, or neither.
    # count(ba.target_name) rather than count(*) - a LEFT JOIN with no match still produces one
    # null-padded row, which count(*) would wrongly count as 1 instead of 0.
    route_coverage = con.execute("""
        SELECT
            d."Substance Name" AS dod_substance_name,
            count(CASE WHEN ba.match_route = 'chembl_inchikey' THEN 1 END) > 0 AS found_via_chembl,
            count(CASE WHEN ba.match_route = 'ctd_inchikey' THEN 1 END) > 0 AS found_via_ctd,
            count(ba.target_name) AS n_rows
        FROM dod.spl d
        LEFT JOIN dod.binding_affinity ba ON ba.dod_substance_name = d."Substance Name"
        GROUP BY d."Substance Name"
        ORDER BY n_rows DESC
    """).df()

    print(f"===== BindingDB Binding Data Per DoD Chemical (chembl_inchikey + ctd_inchikey union) =====")
    print(f"Total BindingDB rows found: {total_rows:,}")
    print(f"Chemicals with data: {(route_coverage['n_rows'] > 0).sum()}/{len(route_coverage)}")
    print(route_coverage)

    neither = route_coverage[~route_coverage['found_via_chembl'] & ~route_coverage['found_via_ctd']]
    print(f"\n===== Found via NEITHER route ({len(neither)}) - candidates for manual lookup =====")
    print(neither['dod_substance_name'].tolist())

    sample = con.execute(f"""
        SELECT dod_substance_name, match_route, target_name, ki_nm, ic50_nm, kd_nm, ec50_nm
        FROM dod.binding_affinity ORDER BY dod_substance_name LIMIT {display_n}
    """).df()
    print(f"\nSample of matched affinity records (first {display_n}):")
    print(sample)

    # Gene -> UniProt accession, human only (same crosswalk confirmed earlier in this project)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _gene_to_uniprot AS
        SELECT DISTINCT csyn.component_synonym AS GeneSymbol, cs.accession AS UniProtAccession
        FROM chembl.component_synonyms csyn
        JOIN chembl.component_sequences cs ON cs.component_id = csyn.component_id
        WHERE csyn.syn_type IN ('GENE_SYMBOL', 'GENE_SYMBOL_OTHER') AND cs.organism ILIKE '%homo sapiens%'
    """)
    overlap = con.execute("""
        SELECT DISTINCT ba.dod_substance_name, ba.target_name, gu.GeneSymbol
        FROM dod.binding_affinity ba
        JOIN _gene_to_uniprot gu ON gu.UniProtAccession IN (ba.uniprot_swissprot, ba.uniprot_trembl)
        JOIN _dod_chem_gene_split dgs ON dgs.GeneSymbol = gu.GeneSymbol AND dgs.Substance = ba.dod_substance_name
    """).df()
    print(f"\n===== Overlap: Same Chemical, Same Gene, in Both CTD and BindingDB =====")
    print(f"{len(overlap)} BindingDB record(s) land on a gene CTD also reports an interaction with for that chemical")
    if len(overlap):
        print(overlap)

    return {"total_rows": total_rows, "route_coverage": route_coverage, "overlap": overlap}


# =========================================================================================
# STEP 4: Per-chemical gene-disease NetworkX graphs
# =========================================================================================
# A simpler, single-chemical sibling of analyze.py's build_type_disease_network - no
# interaction-type filtering across the whole database, no already-inferred/novel tiering
# (that machinery existed there to rank/filter tens of thousands of candidate pairs down to a
# drawable number for one interaction type; here there's only ever one chemical's own genes).
#
# Edges:
#   chemical -> gene    any CTD chem_gene_ixns interaction for this chemical (uniform style -
#                       see note below on why the increases/decreases/affects direction isn't
#                       colored anymore)
#   gene -> disease     direct evidence only, from ctd.genes_diseases
#   chemical -> disease direct evidence only, from ctd.chemicals_diseases
# A disease reached both ways (through a drawn gene AND directly from the chemical) is visually
# apparent as two edges converging on the same node - that convergence is the "gene plausibly
# explains this disease" story, without needing an explicit corroborated/novel color scheme.
#
# Chemical-gene edges are NOT colored by interaction direction (increases/decreases/affects)
# anymore - a gene commonly has both "increases" and "decreases" records for the same chemical
# via different mechanisms, so picking one to color by (the previous code used any_value(),
# an arbitrary pick) doesn't represent anything real. Removed rather than fixed with a "properly
# aggregated" direction, since with only one chemical in view (no interaction-type filtering
# happening here) a single edge color wouldn't add real information either way.
#
# `gene_limit` keeps only the chemical's `gene_limit` most-studied genes - ranked by DISTINCT
# PubMed IDs tying that gene to this chemical in ctd.chem_gene_ixns (breadth of independent
# evidence), not by raw interaction-action row count. The two can diverge a lot: one paper can
# report five distinct actions for a gene in a single publication, which would inflate a
# row-count ranking without reflecting any independent replication. PubMed count is more
# resistant to that (though not immune - a single broad-panel study can still spread PMIDs
# thinly across many genes; no CTD-derived metric fully escapes "what got studied" bias).
# Distinct interaction-action count is still computed and printed alongside for transparency,
# it's just not what genes get ranked/filtered by.
# `diseases_per_gene` caps disease nodes per gene, preferring diseases the chemical is ALSO
# independently linked to (corroborating) before filling remaining slots with the gene's others.

def chemical_gene_disease_network(con, chemical_name, human_only=True, gene_limit=20, diseases_per_gene=2):
    """Builds and returns the chemical -> gene -> disease graph for one chemical (does not
    save/show the figure - see save_all_dod_networks for the per-chemical driver)."""
    organism_filter = "AND Organism ILIKE '%Homo sapiens%'" if human_only else ""

    # Two independent aggregations (not one query unnesting both PubMedIDs and InteractionActions
    # together - that would cross-join the two pipe-delimited lists per row and give nonsense
    # counts), merged afterward on GeneSymbol.
    n_studies = con.execute(f"""
        SELECT GeneSymbol, count(DISTINCT PubMedID) AS n_studies
        FROM (
            SELECT GeneSymbol, UNNEST(string_split(PubMedIDs, '|')) AS PubMedID
            FROM ctd.chem_gene_ixns
            WHERE ChemicalName = ? AND PubMedIDs IS NOT NULL {organism_filter}
        )
        GROUP BY GeneSymbol
    """, [chemical_name]).df()
    n_interactions = con.execute(f"""
        SELECT GeneSymbol, count(DISTINCT InteractionAction) AS n_interactions
        FROM (
            SELECT GeneSymbol, UNNEST(string_split(InteractionActions, '|')) AS InteractionAction
            FROM ctd.chem_gene_ixns
            WHERE ChemicalName = ? AND InteractionActions IS NOT NULL {organism_filter}
        )
        GROUP BY GeneSymbol
    """, [chemical_name]).df()
    gene_evidence = n_studies.merge(n_interactions, on='GeneSymbol', how='outer').fillna(0)
    gene_evidence[['n_studies', 'n_interactions']] = gene_evidence[['n_studies', 'n_interactions']].astype(int)

    # "Death" excluded - see common_diseases for why (non-specific outcome, not a clinical
    # diagnosis, dominates every chemical's disease list without saying anything about mechanism).
    chem_diseases = set(con.execute("""
        SELECT DISTINCT DiseaseName FROM ctd.chemicals_diseases
        WHERE ChemicalName = ? AND DirectEvidence IS NOT NULL AND DiseaseName != 'Death'
    """, [chemical_name]).df()['DiseaseName'])

    if gene_evidence.empty:
        print(f"{chemical_name}: no CTD gene interactions found - skipping")
        return None

    # Keep the gene_limit genes with the most distinct studies (see comment block above).
    gene_evidence = gene_evidence.sort_values('n_studies', ascending=False)
    kept_evidence = gene_evidence.head(gene_limit)
    kept_genes = set(kept_evidence['GeneSymbol'])
    chem_genes = kept_evidence[['GeneSymbol']].reset_index(drop=True)

    print(f"{chemical_name}: top genes by distinct studies (of {len(gene_evidence)} total genes):")
    print(kept_evidence.head(10).to_string(index=False))

    con.register('_chem_genes_df', chem_genes)
    gene_disease = con.execute("""
        SELECT cg.GeneSymbol, gd.DiseaseName
        FROM _chem_genes_df cg
        JOIN ctd.genes_diseases gd ON gd.GeneSymbol = cg.GeneSymbol AND gd.DirectEvidence IS NOT NULL
        WHERE gd.DiseaseName != 'Death'
    """).df()
    con.unregister('_chem_genes_df')

    gene_disease['corroborating'] = gene_disease['DiseaseName'].isin(chem_diseases)
    gene_disease = gene_disease[gene_disease['GeneSymbol'].isin(kept_genes)]

    # Per gene, corroborating diseases first, then fill up to diseases_per_gene with the rest.
    gene_disease = gene_disease.sort_values('corroborating', ascending=False)
    gene_disease = gene_disease.groupby('GeneSymbol', group_keys=False).head(diseases_per_gene)

    print(f"{chemical_name}: {len(kept_genes)} genes (of {len(gene_evidence)} total), "
          f"{gene_disease['DiseaseName'].nunique()} disease nodes, "
          f"{gene_disease['corroborating'].sum()} gene-corroborated disease edges")

    G = nx.Graph()
    node_kind = {chemical_name: 'chemical'}

    for gene in chem_genes['GeneSymbol']:
        node_kind[gene] = 'gene'
        G.add_edge(chemical_name, gene, color="#55A868", style='solid', width=1.5)

    for _, row in gene_disease.iterrows():
        node_kind[row['DiseaseName']] = 'disease'
        G.add_edge(row['GeneSymbol'], row['DiseaseName'], color="#1f77b4", style='dashed', width=1.0)

    for disease in chem_diseases & set(gene_disease['DiseaseName']):
        G.add_edge(chemical_name, disease, color="#ff7f0e", style='solid', width=2.5)

    kind_color = {'chemical': '#4C72B0', 'gene': '#55A868', 'disease': '#C44E52'}
    pos = nx.spring_layout(G, seed=42, k=0.6)
    plt.figure(figsize=(14, 10))

    for kind, color in kind_color.items():
        nodes = [n for n in G.nodes() if node_kind[n] == kind]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=color, node_size=400, label=kind.capitalize())

    for style in ('solid', 'dashed'):
        edges = [(u, v) for u, v in G.edges() if G[u][v]['style'] == style]
        if edges:
            nx.draw_networkx_edges(G, pos, edgelist=edges, style=style,
                                    edge_color=[G[u][v]['color'] for u, v in edges],
                                    width=[G[u][v]['width'] for u, v in edges])

    nx.draw_networkx_labels(G, pos, font_size=7)
    plt.title(f"{chemical_name}: chemical -> gene -> disease network\n"
              f"orange = chemical is ALSO independently linked to that disease (gene-corroborated)")

    node_handles, _ = plt.gca().get_legend_handles_labels()
    edge_handles = [
        Line2D([0], [0], color="#55A868", lw=1.5, label='Chemical-gene interaction'),
        Line2D([0], [0], color="#1f77b4", lw=1.0, linestyle='dashed', label='Gene-disease link'),
        Line2D([0], [0], color="#ff7f0e", lw=2.5, label='Gene-corroborated chemical-disease link'),
    ]
    plt.legend(handles=node_handles + edge_handles, scatterpoints=1, loc='best', fontsize=8)
    plt.axis('off')
    plt.tight_layout()

    return {"graph": G, "gene_evidence": kept_evidence}


def save_all_dod_networks(con, out_dir=None, human_only=True, gene_limit=20, diseases_per_gene=2):
    """Builds and saves one PNG per DoD chemical that has CTD gene interactions - 35 interactive
    windows isn't practical, so figures are written to disk instead of shown. Also persists each
    chemical's top-gene evidence (n_studies/n_interactions) to dod.gene_evidence, so the report
    can show that table without re-deriving it from the image."""
    out_dir = Path(out_dir) if out_dir else NETWORKS_OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    substances = con.execute('SELECT "Substance Name" FROM dod.spl ORDER BY "Substance Name"').df()['Substance Name']

    print(f"===== Building Per-Chemical Gene-Disease Networks =====")
    saved = []
    all_gene_evidence = []
    for substance in substances:
        result = chemical_gene_disease_network(con, substance, human_only=human_only,
                                                gene_limit=gene_limit, diseases_per_gene=diseases_per_gene)
        if result is None:
            plt.close('all')
            continue
        G = result["graph"]
        gene_evidence = result["gene_evidence"].copy()
        gene_evidence.insert(0, 'dod_substance_name', substance)
        all_gene_evidence.append(gene_evidence)

        safe_name = re.sub(r'[^A-Za-z0-9]+', '_', substance).strip('_')
        out_path = out_dir / f"{safe_name}.png"
        plt.savefig(out_path, dpi=150)
        plt.close('all')
        saved.append(out_path)

    combined_evidence = pd.concat(all_gene_evidence, ignore_index=True) if all_gene_evidence else pd.DataFrame()
    if len(combined_evidence):
        con.register('_combined_gene_evidence_df', combined_evidence)
        con.execute("CREATE OR REPLACE TABLE dod.gene_evidence AS SELECT * FROM _combined_gene_evidence_df")
        con.unregister('_combined_gene_evidence_df')

    print(f"\nSaved {len(saved)}/{len(substances)} chemical networks to {out_dir}")
    return saved


# =========================================================================================
# STEP 5: FDA-approved candidate drugs with an opposing mechanism, per DoD chemical
# =========================================================================================
# For every gene a chemical is reported (in CTD) to increase or decrease - ALL of them, not
# just the gene_limit-capped subset drawn in step 4's networks, since this is a cheap targeted
# lookup rather than something that needs capping for a drawing - this looks for an APPROVED
# ChEMBL drug (max_phase=4, not withdrawn) whose curated mechanism-of-action on that same gene's
# protein target runs the OPPOSITE direction: if the chemical increases a gene's activity/
# expression, an antagonist/inhibitor/blocker is a candidate to counteract it; if it decreases
# the gene, an agonist/activator is a candidate. "Opposite" is read directly off ChEMBL's own
# action_type.parent_type classification (POSITIVE MODULATOR vs NEGATIVE MODULATOR) - not
# something hand-rolled here. "affects" interactions are skipped (no clear opposite for a
# non-directional claim), and only genes with a human gene-symbol -> UniProt -> ChEMBL-target
# crosswalk (same crosswalk used throughout this file) can be checked at all.
#
# This is hypothesis generation, not a treatment recommendation. "Same target, opposite
# action" doesn't account for dose, tissue distribution, whether the target is druggable at the
# exposure level involved, or whether the toxic mechanism is even mediated through that gene
# rather than a downstream/secondary effect. Any candidate here needs pharmacology/toxicology
# expert review before being treated as actionable.
#
# Requires match_dod_chemicals(con) to have already run this session (dod.spl).

def opposing_drug_candidates(con, chemical_name, human_only=True, display_n=20):
    """For one chemical: every gene's DOMINANT direction from CTD (increases/decreases, chosen
    by which has more distinct supporting studies - see below), matched to an approved ChEMBL
    drug with the opposite mechanism-of-action on that gene's target. Prints the match rate
    (pairs with a candidate / total directional pairs checked) and returns the candidates."""
    organism_filter = "AND cg.Organism ILIKE '%Homo sapiens%'" if human_only else ""

    # Distinct-study count per (gene, direction), used to pick ONE dominant direction per gene
    # rather than testing every direction CTD ever recorded equally. Confirmed necessary on real
    # data: Sarin/ACHE has 87 studies for "decreases" (the correct, textbook mechanism) vs only 8
    # for "increases" - without this, the minority "increases" claim gets tested just as readily
    # and can surface backwards-looking candidates (e.g. recommending MORE cholinesterase
    # inhibition to "counteract" a cholinesterase-inhibiting nerve agent).
    #
    # Uses a FROM-clause (lateral) UNNEST of both InteractionActions and PubMedIDs so they cross
    # join properly (every action in a row is supported by every citation listed for that same
    # row) - confirmed empirically that unnesting both in the SELECT list instead does a
    # positional zip with NULL-padding on length mismatches, which would silently miscount.
    direction_studies = con.execute(f"""
        SELECT GeneSymbol, lower(trim(split_part(RawAction, '^', 1))) AS Direction,
               count(DISTINCT PubMedID) AS n_studies
        FROM ctd.chem_gene_ixns cg,
             UNNEST(string_split(cg.InteractionActions, '|')) AS a1(RawAction),
             UNNEST(string_split(cg.PubMedIDs, '|')) AS a2(PubMedID)
        WHERE cg.ChemicalName = ? AND cg.InteractionActions IS NOT NULL AND cg.PubMedIDs IS NOT NULL
              {organism_filter}
        GROUP BY GeneSymbol, Direction
    """, [chemical_name]).df()
    direction_studies = direction_studies[direction_studies['Direction'].isin(['increases', 'decreases'])]

    # Keep only the top-studied direction per gene (ties keep both - no principled tiebreak).
    # n_studies is carried forward into the candidates below, used as a secondary ranking key
    # in top_counteracting_drugs - confirmed necessary: many drugs tie at "opposes 1 gene," and
    # without a tiebreak the choice among them is arbitrary (lost Pralidoxime/ACHE for Sarin to
    # an arbitrary pick among 7 equally-"1 gene" candidates on the first pass).
    direction_studies = direction_studies.sort_values('n_studies', ascending=False)
    max_studies = direction_studies.groupby('GeneSymbol')['n_studies'].transform('max')
    directions = direction_studies[direction_studies['n_studies'] == max_studies][['GeneSymbol', 'Direction', 'n_studies']]

    if directions.empty:
        print(f"{chemical_name}: no directional (increases/decreases) gene interactions found - skipping")
        return {"chemical_name": chemical_name, "n_pairs_checked": 0, "n_pairs_with_candidate": 0,
                "candidates": pd.DataFrame()}

    # Gene symbol -> ChEMBL target (tid), human only - same crosswalk used throughout this file,
    # extended one hop further from UniProt accession to the target_dictionary row it belongs to.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _gene_to_tid AS
        SELECT DISTINCT csyn.component_synonym AS GeneSymbol, td.tid,
               td.chembl_id AS target_chembl_id, td.pref_name AS target_pref_name
        FROM chembl.component_synonyms csyn
        JOIN chembl.component_sequences cs ON cs.component_id = csyn.component_id
        JOIN chembl.target_components tc ON tc.component_id = cs.component_id
        JOIN chembl.target_dictionary td ON td.tid = tc.tid
        WHERE csyn.syn_type IN ('GENE_SYMBOL', 'GENE_SYMBOL_OTHER') AND cs.organism ILIKE '%homo sapiens%'
    """)

    opposite_parent_type = {'increases': 'NEGATIVE MODULATOR', 'decreases': 'POSITIVE MODULATOR'}
    all_candidates = []
    for direction, wanted_parent_type in opposite_parent_type.items():
        genes = directions.loc[directions['Direction'] == direction, ['GeneSymbol', 'n_studies']]
        if genes.empty:
            continue
        con.register('_genes_df', genes)
        candidates = con.execute("""
            SELECT gt.GeneSymbol, g.n_studies, gt.target_chembl_id, gt.target_pref_name,
                   md.chembl_id AS drug_chembl_id, md.pref_name AS drug_name,
                   md.withdrawn_flag, md.black_box_warning,
                   dm.action_type, dm.mechanism_of_action
            FROM _genes_df g
            JOIN _gene_to_tid gt ON gt.GeneSymbol = g.GeneSymbol
            JOIN chembl.drug_mechanism dm ON dm.tid = gt.tid
            JOIN chembl.action_type atype ON atype.action_type = dm.action_type
            JOIN chembl.molecule_dictionary md ON md.molregno = dm.molregno
            WHERE atype.parent_type = ? AND md.max_phase = 4 AND md.withdrawn_flag = 0
        """, [wanted_parent_type]).df()
        con.unregister('_genes_df')
        if not candidates.empty:
            candidates.insert(0, 'ctd_direction', direction)
            candidates.insert(1, 'wanted_drug_action', wanted_parent_type)
            all_candidates.append(candidates)

    all_candidates = pd.concat(all_candidates, ignore_index=True) if all_candidates else pd.DataFrame()

    n_pairs_checked = len(directions)
    n_pairs_with_candidate = (
        all_candidates[['GeneSymbol', 'ctd_direction']].drop_duplicates().shape[0] if len(all_candidates) else 0
    )
    print(f"===== {chemical_name}: Opposing-Mechanism Drug Candidates =====")
    print(f"{n_pairs_with_candidate}/{n_pairs_checked} (gene, direction) pairs have an approved-drug "
          f"candidate with an opposing mechanism")
    if len(all_candidates):
        print(all_candidates.head(display_n))

    return {"chemical_name": chemical_name, "n_pairs_checked": n_pairs_checked,
            "n_pairs_with_candidate": n_pairs_with_candidate, "candidates": all_candidates}


def find_all_opposing_drug_candidates(con, human_only=True):
    """Runs opposing_drug_candidates for every DoD chemical, writes the combined result to
    dod.opposing_drug_candidates, and prints the overall match rate across the whole list."""
    substances = con.execute('SELECT "Substance Name" FROM dod.spl ORDER BY "Substance Name"').df()['Substance Name']

    all_results = []
    total_checked = 0
    total_with_candidate = 0
    for substance in substances:
        result = opposing_drug_candidates(con, substance, human_only=human_only)
        total_checked += result["n_pairs_checked"]
        total_with_candidate += result["n_pairs_with_candidate"]
        if len(result["candidates"]):
            tagged = result["candidates"].copy()
            tagged.insert(0, 'dod_substance_name', substance)
            all_results.append(tagged)

    combined = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    if len(combined):
        con.register('_combined_candidates_df', combined)
        con.execute("CREATE OR REPLACE TABLE dod.opposing_drug_candidates AS SELECT * FROM _combined_candidates_df")
        con.unregister('_combined_candidates_df')

    print(f"\n===== Overall: Opposing-Mechanism Drug Candidates Across the DoD List =====")
    print(f"{total_with_candidate}/{total_checked} (chemical, gene, direction) pairs have an approved-drug "
          f"candidate with an opposing mechanism ({total_with_candidate / total_checked * 100:.1f}%)"
          if total_checked else "No directional gene interactions found across the DoD list")
    print(f"Distinct chemicals with at least one candidate: {combined['dod_substance_name'].nunique() if len(combined) else 0}")

    return combined


def top_counteracting_drugs(con, chemical_name, top_n=3):
    """The actual per-chemical recommendation list: aggregates dod.opposing_drug_candidates
    BY DRUG (rather than by gene) for one chemical, ranking drugs by how many of the chemical's
    distinct CTD-affected genes each drug opposes. A drug opposing 5 of a chemical's genes is a
    more broadly-relevant candidate than one that only opposes 1 - "most common opposite effect
    across everything this chemical affects", not just a flat per-gene match list.

    NEGATIVE MODULATOR and POSITIVE MODULATOR candidates are ranked SEPARATELY (top_n of each,
    not one combined top_n) rather than pooled into a single ranking. A combined ranking lets
    whichever direction happens to have broader multi-gene coverage crowd out the other
    entirely - confirmed concretely with Sarin: Pralidoxime (the real antidote) only opposes
    ACHE (1 gene) and got pushed off a combined top-5 list dominated by unrelated drugs that
    each opposed 2 genes. Splitting guarantees both directions get a chance to appear.

    Ties on n_genes_countered are broken by max(n_studies) - the strongest per-gene evidence
    (distinct PubMed count, carried through from opposing_drug_candidates) among the genes a
    drug opposes. Also confirmed necessary: for Sarin, 7 different POSITIVE MODULATOR drugs
    each oppose exactly 1 gene, so an untiebroken top-3 picked an arbitrary 3 of the 7 and
    initially missed Pralidoxime (ACHE, 37 supporting studies) in favor of candidates whose
    genes had 1 or zero supporting studies.

    Requires find_all_opposing_drug_candidates(con) to have already run this session
    (dod.opposing_drug_candidates)."""
    ranked = con.execute("""
        SELECT * FROM (
            SELECT wanted_drug_action, drug_chembl_id, drug_name,
                   count(DISTINCT GeneSymbol) AS n_genes_countered,
                   array_agg(DISTINCT GeneSymbol) AS genes_countered,
                   max(n_studies) AS max_gene_n_studies,
                   bool_or(black_box_warning = 1) AS has_black_box_warning,
                   row_number() OVER (
                       PARTITION BY wanted_drug_action
                       ORDER BY count(DISTINCT GeneSymbol) DESC, max(n_studies) DESC
                   ) AS rnk
            FROM dod.opposing_drug_candidates
            WHERE dod_substance_name = ?
            GROUP BY wanted_drug_action, drug_chembl_id, drug_name
        )
        WHERE rnk <= ?
        ORDER BY wanted_drug_action, rnk
    """, [chemical_name, top_n]).df()
    if len(ranked):
        ranked = ranked.drop(columns=['rnk'])

    print(f"===== {chemical_name}: Top {top_n} Counteracting Drug Candidates (per direction) =====")
    if ranked.empty:
        print("(no approved-drug candidates found for this chemical)")
    else:
        print(ranked)

    return ranked


def find_all_top_counteracting_drugs(con, top_n=3):
    """Builds the top_counteracting_drugs list for every DoD chemical and writes the combined
    result to dod.top_counteracting_drugs - the actual per-chemical recommendation table."""
    substances = con.execute('SELECT "Substance Name" FROM dod.spl ORDER BY "Substance Name"').df()['Substance Name']

    all_results = []
    for substance in substances:
        ranked = top_counteracting_drugs(con, substance, top_n=top_n)
        if not ranked.empty:
            ranked = ranked.copy()
            ranked.insert(0, 'dod_substance_name', substance)
            all_results.append(ranked)

    combined = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    if len(combined):
        con.register('_combined_top_drugs_df', combined)
        con.execute("CREATE OR REPLACE TABLE dod.top_counteracting_drugs AS SELECT * FROM _combined_top_drugs_df")
        con.unregister('_combined_top_drugs_df')

    print(f"\n===== Overall: Top Counteracting Drugs Table =====")
    print(f"{combined['dod_substance_name'].nunique() if len(combined) else 0}/{len(substances)} chemicals "
          f"have at least one recommended counteracting drug")

    return combined


# =========================================================================================
# STEP 6: Aggregated HTML report - one summary section + one section per chemical
# =========================================================================================
# A single document combining everything the earlier steps computed: overview stats and
# commonalities (steps 2), then per chemical: ChEMBL match/mechanism info (step 1), BindingDB
# binding data (step 3), top genes by evidence (step 4), the network image (step 4), the
# chemical's own CTD-linked diseases, and the top counteracting drug candidates (step 5).
#
# "Death" is excluded from every disease listing here too (see common_diseases/
# chemical_gene_disease_network for why - it's CTD's most common "shared disease" for almost
# any toxicant but carries no clinical/mechanistic specificity).
#
# Requires steps 1-5 to have already run this session (reads dod.chembl_match,
# dod.chembl_activities, dod.binding_affinity, dod.gene_evidence, dod.top_counteracting_drugs,
# plus the saved PNGs in data/dod/networks/). An HTML report (not plain text, unlike
# analyze.py/report.py's ATSDR sibling) since this one needs to actually display the network
# images, not just point at their file paths.

REPORT_OUT_PATH = PROJECT_ROOT / 'data' / 'dod' / 'dod_report.html'

REPORT_CSS = """
<style>
  body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; max-width: 1100px;
         margin: 2em auto; padding: 0 1em; line-height: 1.5; color: #1a1a1a; }
  h1 { border-bottom: 3px solid #4C72B0; padding-bottom: 0.3em; }
  h2 { margin-top: 3em; border-bottom: 1px solid #ccc; padding-bottom: 0.2em; color: #4C72B0; }
  h3 { margin-top: 1.5em; color: #333; }
  table { border-collapse: collapse; margin: 0.8em 0; font-size: 0.9em; }
  th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
  th { background: #f0f2f5; }
  img { max-width: 100%; height: auto; border: 1px solid #ccc; margin: 1em 0; }
  .meta { color: #666; font-size: 0.9em; }
  .empty { color: #888; font-style: italic; }
  .chem-section { margin-bottom: 4em; }

  /* Standard portrait page. Tables are kept narrow (few, consolidated columns - see
     chembl_activities/binding queries below) so they fit without landscape orientation.
     table-layout: auto (not fixed) lets column widths follow content instead of forcing every
     column to the same width, which is what made tables look cramped/"square" before; width:
     100% + word-wrap still caps the table at the page width so nothing overflows. */
  @media print {
    @page { size: portrait; margin: 1.2cm; }
    body { max-width: 100%; font-size: 0.85em; }
    table { width: 100%; table-layout: auto; font-size: 0.8em; }
    th, td { word-wrap: break-word; overflow-wrap: break-word; white-space: normal; }
    img { max-width: 90%; }
  }
</style>
"""


def _df_html(df, empty_msg="(no data found)", max_rows=25):
    if df is None or len(df) == 0:
        return f'<p class="empty">{empty_msg}</p>'
    note = ""
    if len(df) > max_rows:
        note = f'<p class="meta">Showing {max_rows} of {len(df)} rows.</p>'
        df = df.head(max_rows)
    # na_rep alone doesn't reliably blank Python None in DuckDB-sourced object columns -
    # replace explicitly rather than rely on to_html's own NA detection.
    df = df.where(df.notna(), "")
    return note + df.to_html(index=False, na_rep="")


def generate_dod_report(con, out_path=None):
    out_path = Path(out_path) if out_path else REPORT_OUT_PATH
    html = ['<!doctype html><html><head><meta charset="utf-8"><title>DoD Chemicals of Concern Report</title>', REPORT_CSS, '</head><body>']
    html.append("<h1>DoD/VA Chemicals of Concern - Analysis Report</h1>")

    substances = con.execute('SELECT "Substance Name", "CASRN" FROM dod.spl ORDER BY "Substance Name"').df()

    # ---------------- Summary section ----------------
    html.append("<h2>Summary</h2>")

    n_total = len(substances)
    n_chembl_matched = con.execute("SELECT count(*) FROM dod.chembl_match WHERE matched").fetchone()[0]
    n_with_activity = con.execute("SELECT count(DISTINCT dod_substance_name) FROM dod.chembl_activities").fetchone()[0]
    n_with_binding = con.execute("SELECT count(DISTINCT dod_substance_name) FROM dod.binding_affinity").fetchone()[0]
    n_with_drug_candidate = con.execute("SELECT count(DISTINCT dod_substance_name) FROM dod.top_counteracting_drugs").fetchone()[0]
    html.append(f"""
        <ul>
          <li>{n_total} chemicals in the DoD list</li>
          <li>{n_chembl_matched}/{n_total} matched to a ChEMBL compound</li>
          <li>{n_with_activity}/{n_total} have ChEMBL bioactivity/mechanism data</li>
          <li>{n_with_binding}/{n_total} have BindingDB binding-affinity data (chembl_inchikey + ctd_inchikey routes combined)</li>
          <li>{n_with_drug_candidate}/{n_total} have at least one approved-drug counteracting candidate</li>
        </ul>
    """)

    html.append("<h3>Genes shared across the DoD list (hub genes)</h3>")
    hub = common_genes(con)
    hub_df = hub["hub_genes"].reset_index()
    hub_df.columns = ["GeneSymbol", "n_chemicals"]
    html.append(_df_html(hub_df, max_rows=20))

    html.append("<h3>Diseases shared across the DoD list</h3>")
    diseases = common_diseases(con)
    disease_df = diseases["shared_diseases"].reset_index()
    disease_df.columns = ["DiseaseName", "n_chemicals"]
    html.append(_df_html(disease_df, max_rows=20))

    # ---------------- Per-chemical sections ----------------
    for _, row in substances.iterrows():
        substance, casrn = row["Substance Name"], row["CASRN"]
        html.append(f'<div class="chem-section">')
        html.append(f"<h2>{substance}</h2>")
        if casrn:
            html.append(f'<p class="meta">CASRN: {casrn}</p>')

        # ChEMBL match status
        match_info = con.execute("""
            SELECT cm.matched, cm.match_method, md.chembl_id, md.pref_name AS chembl_pref_name
            FROM dod.chembl_match cm
            LEFT JOIN chembl.molecule_dictionary md ON md.molregno = cm.chembl_molregno
            WHERE cm.substance_name = ?
        """, [substance]).df()
        html.append("<h3>ChEMBL Match</h3>")
        if len(match_info) and match_info.iloc[0]["matched"]:
            m = match_info.iloc[0]
            pref_name_suffix = f" ({m['chembl_pref_name']})" if m['chembl_pref_name'] else ""
            html.append(f"<p>Matched via <b>{m['match_method']}</b> to {m['chembl_id']}{pref_name_suffix}</p>")
        else:
            html.append('<p class="empty">No ChEMBL match found for this substance.</p>')

        # ChEMBL bioactivity table dropped (was here) - checked the real data: 0/1274 rows had
        # mechanism_of_action populated for this DoD-specific chemical set, 36% were blank noise
        # (assay run, no value, no mechanism - e.g. "ADMET/Bacterial Biotransformation"), and most
        # of the rest were non-mechanistic (LogP, rat toxicology panel readouts, mortality %).
        # The genuinely useful chemical-target relationships are already covered by the network
        # diagram (gene evidence) and the BindingDB affinity table below.

        # BindingDB binding affinity. Ki/IC50/Kd/EC50 are combined into one "affinity" column -
        # confirmed only one of the four is ever populated per row (295/299 rows have exactly
        # one; the rest have none), so nothing is lost by not giving each its own column.
        binding = con.execute("""
            SELECT target_name, match_route,
                   CASE WHEN ki_nm IS NOT NULL THEN 'Ki: ' || ki_nm || ' nM'
                        WHEN ic50_nm IS NOT NULL THEN 'IC50: ' || ic50_nm || ' nM'
                        WHEN kd_nm IS NOT NULL THEN 'Kd: ' || kd_nm || ' nM'
                        WHEN ec50_nm IS NOT NULL THEN 'EC50: ' || ec50_nm || ' nM'
                        ELSE NULL END AS affinity,
                   source
            FROM dod.binding_affinity WHERE dod_substance_name = ?
            ORDER BY target_name
        """, [substance]).df()
        html.append("<h3>BindingDB Binding Affinity</h3>")
        html.append(_df_html(binding, empty_msg="No BindingDB binding-affinity data found (see Step 3 notes on coverage limitations).", max_rows=15))

        # Top genes by evidence
        gene_evidence = con.execute("""
            SELECT GeneSymbol, n_studies, n_interactions
            FROM dod.gene_evidence WHERE dod_substance_name = ?
            ORDER BY n_studies DESC
        """, [substance]).df()
        html.append("<h3>Top CTD Genes (by distinct supporting studies)</h3>")
        html.append(_df_html(gene_evidence, empty_msg="No CTD gene interactions found.", max_rows=15))

        # Network image
        html.append("<h3>Chemical -&gt; Gene -&gt; Disease Network</h3>")
        safe_name = re.sub(r'[^A-Za-z0-9]+', '_', substance).strip('_')
        png_path = NETWORKS_OUT_DIR / f"{safe_name}.png"
        if png_path.exists():
            html.append(f'<img src="networks/{safe_name}.png" alt="{substance} network">')
        else:
            html.append('<p class="empty">No network image available (no CTD gene interactions).</p>')

        # This chemical's own CTD diseases (direct evidence, Death excluded)
        own_diseases = con.execute("""
            SELECT DISTINCT DiseaseName FROM ctd.chemicals_diseases
            WHERE ChemicalName = ? AND DirectEvidence IS NOT NULL AND DiseaseName != 'Death'
            ORDER BY DiseaseName
        """, [substance]).df()
        html.append("<h3>Diseases Linked to This Chemical (direct evidence, excl. Death)</h3>")
        html.append(_df_html(own_diseases, empty_msg="No direct-evidence disease links found."))

        # Top counteracting drug candidates - negative/positive modulators ranked separately
        # (see top_counteracting_drugs docstring for why: a combined ranking can let one
        # direction's broader multi-gene coverage crowd out the other entirely).
        top_drugs = con.execute("""
            SELECT wanted_drug_action, drug_name, drug_chembl_id, n_genes_countered, genes_countered,
                   max_gene_n_studies, has_black_box_warning
            FROM dod.top_counteracting_drugs WHERE dod_substance_name = ?
            ORDER BY wanted_drug_action, n_genes_countered DESC, max_gene_n_studies DESC
        """, [substance]).df()
        if len(top_drugs):
            top_drugs['genes_countered'] = top_drugs['genes_countered'].apply(lambda a: ', '.join(a))
        html.append("<h3>Top Counteracting Drug Candidates (approved, opposing mechanism, 3 per direction)</h3>")
        html.append(_df_html(top_drugs, empty_msg="No approved-drug counteracting candidates found."))
        html.append("</div>")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    html.append("</body></html>")
    out_path.write_text("\n".join(html), encoding='utf-8')
    print(f"Report written to {out_path}")
    return out_path


if __name__ == '__main__':
    con = connect()

    match_dod_chemicals(con)
    print()
    common_genes(con)
    print()
    common_diseases(con)
    print()
    binding_affinity_per_chemical(con)
    print()
    save_all_dod_networks(con)
    print()
    find_all_opposing_drug_candidates(con)
    print()
    find_all_top_counteracting_drugs(con)
    print()
    generate_dod_report(con)

    con.close()
