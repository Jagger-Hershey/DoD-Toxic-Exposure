# parse_ctd.py
# Author: Jagger Hershey
#
# Script for ingesting bulk gzipped csv files from the Comparative Toxicogenomics Database (CTD)
# directly into the unified DuckDB database. DuckDB reads gzipped, comment-prefixed CSVs
# out-of-core, so this replaces the previous pandas chunked-read + parquet-write pipeline
# (which pandas couldn't even load back afterwards for the 125M-row genes_diseases file).
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'common'))
from db import connect, PROJECT_ROOT

RAW_DIR = PROJECT_ROOT / 'data' / 'bulk' / 'raw'

# (table name, source file, column names, column types) - column order/names/types match
# CTD's documented field order; the raw files ship with no header row.
TABLES = [
    (
        "chem_gene_ixns",
        "CTD_chem_gene_ixns.csv.gz",
        ["ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID", "GeneForms",
         "Organism", "OrganismID", "Interaction", "InteractionActions", "PubMedIDs"],
        ["VARCHAR"] * 11,
    ),
    (
        "chemicals_diseases",
        "CTD_chemicals_diseases.csv.gz",
        ["ChemicalName", "ChemicalID", "CasRN", "DiseaseName", "DiseaseID", "DirectEvidence",
         "InferenceGeneSymbol", "InferenceScore", "OmimIDs", "PubMedIDs"],
        ["VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "DOUBLE", "VARCHAR", "VARCHAR"],
    ),
    (
        "genes_pathways",
        "CTD_genes_pathways.csv.gz",
        ["GeneSymbol", "GeneID", "PathwayName", "PathwayID"],
        ["VARCHAR"] * 4,
    ),
    (
        "genes_diseases",
        "CTD_genes_diseases.csv.gz",
        ["GeneSymbol", "GeneID", "DiseaseName", "DiseaseID", "DirectEvidence",
         "InferenceChemicalName", "InferenceScore", "OmimIDs", "PubMedIDs"],
        ["VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "DOUBLE", "VARCHAR", "VARCHAR"],
    ),
]


# CTD's master chemical vocabulary (name/CasRN/InChIKey/etc per chemical). Parsed separately
# from TABLES above because, unlike the interaction/association files, it has quoted fields
# (chemical names and synonym lists can contain literal commas) and needs quote/escape set
# explicitly, plus ignore_errors/null_padding since a handful of rows don't fully comply.
CHEMICALS_FILE = "CTD_chemicals.csv.gz"
CHEMICALS_COLUMNS = ["ChemicalName", "ChemicalID", "CasRN", "PubChemCID", "PubChemSID", "DTXSID",
                     "InChIKey", "Definition", "ParentIDs", "TreeNumbers", "ParentTreeNumbers",
                     "MESHSynonyms", "CTDCuratedSynonyms"]


def parse_bulk():
    con = connect()

    print(f"{'Table':30s} {'Rows':>14s}")
    print("-" * 46)
    for table, filename, names, types in TABLES:
        path = RAW_DIR / filename
        con.execute(f"""
            CREATE OR REPLACE TABLE ctd.{table} AS
            SELECT * FROM read_csv('{path.as_posix()}', header=false, comment='#', delim=',',
                                    names={names!r}, types={types!r})
        """)
        rows = con.execute(f"SELECT count(*) FROM ctd.{table}").fetchone()[0]
        print(f"ctd.{table:26s} {rows:>14,}")

    chemicals_path = RAW_DIR / CHEMICALS_FILE
    types = ["VARCHAR"] * len(CHEMICALS_COLUMNS)
    con.execute(f"""
        CREATE OR REPLACE TABLE ctd.chemicals AS
        SELECT * FROM read_csv('{chemicals_path.as_posix()}', header=false, comment='#', delim=',',
                                quote='"', escape='"', null_padding=true, ignore_errors=true,
                                names={CHEMICALS_COLUMNS!r}, types={types!r})
    """)
    rows = con.execute("SELECT count(*) FROM ctd.chemicals").fetchone()[0]
    print(f"ctd.{'chemicals':26s} {rows:>14,}")

    con.close()


if __name__ == '__main__':
    parse_bulk()
