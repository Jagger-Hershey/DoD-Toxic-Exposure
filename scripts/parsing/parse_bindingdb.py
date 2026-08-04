# parse_bindingdb.py
# Author: Jagger Hershey
#
# Script for ingesting the BindingDB Ligand-Target-Affinity dataset directly into the
# unified DuckDB database.
#
# BindingDB's "All" export repeats an 11-column block ("...Target Chain N") once per protein
# chain in a target, up to N=50, to cover rare multi-subunit crystal structures. Chain 1 is
# populated for every row; chain 2+ is >94% empty. We keep only chain 1 (target identity/sequence)
# plus the shared ligand/affinity/provenance columns, and record the chain count instead of
# carrying the long tail of near-empty columns.
#
# DuckDB can't read inside a .zip archive, so the source tsv is extracted to a temp file first
# and deleted immediately after ingestion (it's ~9GB uncompressed - not worth keeping alongside
# the 600MB zip it came from).
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'common'))
from db import connect, PROJECT_ROOT

RAW_DIR = PROJECT_ROOT / 'data' / 'bindingdb' / 'raw'
ZIP_FILE = RAW_DIR / 'BindingDB_All_202607_tsv.zip'

# Columns kept as-is: ligand identity, affinity measurements, assay/provenance, and chain-1 target identity.
KEPT_COLUMNS = [
    "BindingDB Reactant_set_id",
    "Ligand SMILES",
    "Ligand InChI",
    "Ligand InChI Key",
    "BindingDB MonomerID",
    "BindingDB Ligand Name",
    "Target Name",
    "Target Source Organism According to Curator or DataSource",
    "Ki (nM)",
    "IC50 (nM)",
    "Kd (nM)",
    "EC50 (nM)",
    "kon (M-1-s-1)",
    "koff (s-1)",
    "pH",
    "Temp (C)",
    "Curation/DataSource",
    "Article DOI",
    "BindingDB Entry DOI",
    "PMID",
    "PubChem AID",
    "Patent Number",
    "Link to Ligand in BindingDB",
    "Link to Target in BindingDB",
    "Link to Ligand-Target Pair in BindingDB",
    "PubChem CID",
    "PubChem SID",
    "ChEBI ID of Ligand",
    "ChEMBL ID of Ligand",
    "DrugBank ID of Ligand",
    "KEGG ID of Ligand",
    "ZINC ID of Ligand",
    "Number of Protein Chains in Target (>1 implies a multichain complex)",
    "BindingDB Target Chain Sequence 1",
    "UniProt (SwissProt) Recommended Name of Target Chain 1",
    "UniProt (SwissProt) Primary ID of Target Chain 1",
    "UniProt (TrEMBL) Primary ID of Target Chain 1",
]


def _find_tsv_member(zf: zipfile.ZipFile) -> str:
    tsv_members = [name for name in zf.namelist() if name.lower().endswith('.tsv')]
    if len(tsv_members) != 1:
        raise ValueError(f"Expected exactly one .tsv file inside zip, found: {tsv_members}")
    return tsv_members[0]


def _extract_tsv() -> Path:
    with zipfile.ZipFile(ZIP_FILE) as zf:
        member_name = _find_tsv_member(zf)
        extracted_path = RAW_DIR / Path(member_name).name
        if not extracted_path.exists():
            print(f"Extracting {member_name} ({zf.getinfo(member_name).file_size / 1e9:.1f} GB) ...")
            zf.extract(member_name, path=RAW_DIR)
    return extracted_path


def parse_bindingdb():
    tsv_path = _extract_tsv()
    con = connect()

    quoted_cols = ", ".join(f'"{c}"' for c in KEPT_COLUMNS)
    n_chains_col = '"Number of Protein Chains in Target (>1 implies a multichain complex)"'

    con.execute(f"""
        CREATE OR REPLACE TABLE bindingdb.activities AS
        SELECT {quoted_cols}
        FROM read_csv('{tsv_path.as_posix()}', delim='\t', header=true, quote='', null_padding=true, all_varchar=true)
    """)

    total_rows = con.execute("SELECT count(*) FROM bindingdb.activities").fetchone()[0]
    multichain_rows = con.execute(f"""
        SELECT count(*) FROM bindingdb.activities
        WHERE {n_chains_col} IS NOT NULL AND {n_chains_col} != '1'
    """).fetchone()[0]

    con.close()
    tsv_path.unlink()  # ~9GB extracted file - not worth keeping alongside the 600MB zip

    print(f"Done: bindingdb.activities")
    print(f"Total rows: {total_rows:,}")
    print(f"Columns kept: {len(KEPT_COLUMNS)}")
    print(f"Rows with dropped multichain (chain 2+) data: {multichain_rows:,}")


if __name__ == '__main__':
    parse_bindingdb()
