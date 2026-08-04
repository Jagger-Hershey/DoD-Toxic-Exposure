# parse_chembl.py
# Author: Jagger Hershey
#
# Script for extracting chemical-target bioactivity and mechanism-of-action data from ChEMBL,
# filtered to the ATSDR 2025 Substance Priority List, into the unified DuckDB database.
# ChEMBL itself is never copied in - it's queried directly from its own SQLite file via the
# `chembl` schema attached read-only by scripts/common/db.py.
#
# ChEMBL compounds aren't keyed by CAS number, so each ATSDR substance is matched to a ChEMBL
# molecule in three passes, most confident first:
#   1. CAS number recorded as a molecule synonym (digit-for-digit, ignoring punctuation)
#   2. Exact substance name match (against molecule_dictionary.pref_name or a synonym)
#   3. CAS number -> InChIKey via PubChem, matched against ChEMBL's own standard_inchi_key
# Unmatched substances (common for elements/inorganics, which ChEMBL doesn't cover) are recorded
# in atsdr.chembl_match rather than silently dropped.
import re
import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'common'))
from db import connect, PROJECT_ROOT

ATSDR_CSV = PROJECT_ROOT / 'data' / 'atsdr' / 'ATSDR-2025-Official-SPL.csv'
PUBCHEM_CACHE_PATH = PROJECT_ROOT / 'data' / 'chembl' / 'raw' / 'pubchem_casrn_inchikey_cache.json'

PUBCHEM_CID_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/{casrn}/cids/JSON"
PUBCHEM_INCHIKEY_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON"
PUBCHEM_REQUEST_DELAY_S = 0.25  # stay under PubChem's ~5 requests/second guidance

# Bioactivity + mechanism-of-action data for every matched molecule in one shot, joined across:
# compound structure -> activities -> assay -> target -> target's protein component (UniProt)
# -> curated mechanism-of-action (drug_mechanism), where one exists for that compound/target pair.
ACTIVITY_JOIN_SQL = """
CREATE OR REPLACE TABLE atsdr.chembl_activities AS
SELECT
    mm.atsdr_substance_name,
    mm.atsdr_casrn,
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


def load_atsdr_spl(con):
    con.execute(f"""
        CREATE OR REPLACE TABLE atsdr.spl AS
        SELECT * FROM read_csv('{ATSDR_CSV.as_posix()}', header=true)
    """)
    return con.execute('SELECT "Substance Name", "CASRN" FROM atsdr.spl').fetchall()


def build_cas_and_name_indexes(con, atsdr_list):
    """Single pass over molecule_synonyms/molecule_dictionary (via the attached chembl schema),
    matched against the (small) ATSDR list, instead of one query per substance against
    multi-million-row tables."""
    casdigits_to_casrn = {re.sub(r'\D', '', casrn): casrn for _, casrn in atsdr_list}
    name_to_casrn = {_normalize_name(name): casrn for name, casrn in atsdr_list if name}

    cas_matches = {}
    name_matches = {}

    for molregno, synonym in con.execute("SELECT molregno, synonyms FROM chembl.molecule_synonyms WHERE synonyms IS NOT NULL").fetchall():
        digits = re.sub(r'\D', '', synonym)
        if digits and digits in casdigits_to_casrn:
            cas_matches.setdefault(casdigits_to_casrn[digits], molregno)
        norm = _normalize_name(synonym)
        if norm and norm in name_to_casrn:
            name_matches.setdefault(name_to_casrn[norm], molregno)

    for molregno, pref_name in con.execute("SELECT molregno, pref_name FROM chembl.molecule_dictionary WHERE pref_name IS NOT NULL").fetchall():
        norm = _normalize_name(pref_name)
        if norm and norm in name_to_casrn:
            name_matches.setdefault(name_to_casrn[norm], molregno)

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


def match_remaining_via_pubchem(con, atsdr_list, already_matched: set):
    remaining = [(name, casrn) for name, casrn in atsdr_list if casrn not in already_matched]
    if not remaining:
        return {}

    cache = json.load(open(PUBCHEM_CACHE_PATH)) if PUBCHEM_CACHE_PATH.exists() else {}

    print(f"Resolving {len(remaining)} unmatched substances via PubChem CASRN -> InChIKey ...")
    inchikey_to_casrn = {}
    for _, casrn in remaining:
        inchikey = _pubchem_inchikey_for_casrn(casrn, cache)
        if inchikey:
            inchikey_to_casrn[inchikey] = casrn
    with open(PUBCHEM_CACHE_PATH, 'w') as f:
        json.dump(cache, f, indent=2)

    if not inchikey_to_casrn:
        return {}

    matches = {}
    for molregno, inchikey in con.execute("SELECT molregno, standard_inchi_key FROM chembl.compound_structures WHERE standard_inchi_key IS NOT NULL").fetchall():
        if inchikey in inchikey_to_casrn:
            matches.setdefault(inchikey_to_casrn[inchikey], molregno)
    return matches


def parse_chembl():
    con = connect()

    atsdr_list = [(name.strip(), casrn.strip()) for name, casrn in load_atsdr_spl(con)]
    name_by_casrn = {casrn: name for name, casrn in atsdr_list}

    cas_matches, name_matches = build_cas_and_name_indexes(con, atsdr_list)
    matches = {casrn: (molregno, 'cas_synonym') for casrn, molregno in cas_matches.items()}
    for casrn, molregno in name_matches.items():
        matches.setdefault(casrn, (molregno, 'name_match'))

    pubchem_matches = match_remaining_via_pubchem(con, atsdr_list, set(matches))
    for casrn, molregno in pubchem_matches.items():
        matches.setdefault(casrn, (molregno, 'pubchem_inchikey'))

    # atsdr.chembl_match: every ATSDR substance, matched or not
    con.execute("""
        CREATE OR REPLACE TABLE atsdr.chembl_match (
            substance_name VARCHAR, casrn VARCHAR, chembl_molregno BIGINT,
            match_method VARCHAR, matched BOOLEAN
        )
    """)
    match_rows = [
        (name, casrn, matches[casrn][0] if casrn in matches else None,
         matches[casrn][1] if casrn in matches else None, casrn in matches)
        for name, casrn in atsdr_list
    ]
    con.executemany(
        "INSERT INTO atsdr.chembl_match VALUES (?, ?, ?, ?, ?)", match_rows
    )

    # One SQL join across the attached chembl schema for every matched compound at once,
    # instead of looping per-molregno in Python.
    con.execute("CREATE OR REPLACE TEMP TABLE _matched_molregnos (molregno BIGINT, atsdr_substance_name VARCHAR, atsdr_casrn VARCHAR, match_method VARCHAR)")
    con.executemany(
        "INSERT INTO _matched_molregnos VALUES (?, ?, ?, ?)",
        [(molregno, name_by_casrn[casrn], casrn, method) for casrn, (molregno, method) in matches.items()],
    )
    con.execute(ACTIVITY_JOIN_SQL)

    total_rows = con.execute("SELECT count(*) FROM atsdr.chembl_activities").fetchone()[0]
    substances_with_activity = con.execute("SELECT count(DISTINCT atsdr_casrn) FROM atsdr.chembl_activities").fetchone()[0]
    method_counts = dict(con.execute("SELECT match_method, count(*) FROM atsdr.chembl_match WHERE matched GROUP BY match_method").fetchall())

    con.close()

    print(f"Done: atsdr.chembl_activities, atsdr.chembl_match")
    print(f"Matched {len(matches)}/{len(atsdr_list)} ATSDR substances to a ChEMBL compound")
    print(f"  of which {substances_with_activity} had activity/bioactivity data")
    print(f"Match methods: {method_counts}")
    print(f"Total activity rows written: {total_rows:,}")


if __name__ == '__main__':
    parse_chembl()
