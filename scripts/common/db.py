# db.py
# Author: Jagger Hershey
#
# Shared DuckDB connection for the project's unified database. ChEMBL stays in its own
# 30GB SQLite file (never copied into DuckDB's storage) and is attached read-only on every
# connection so it's directly joinable against the native ctd/bindingdb/atsdr/dod tables.
#
# atsdr.* holds the original ATSDR 2025 Substance Priority List pipeline (unchanged).
# dod.* holds the separate DoD/VA chemicals-of-concern list (see parse_chembl_dod.py) - kept
# as its own schema rather than replacing atsdr, since the ATSDR-specific ranked analysis in
# analyze.py/report.py still relies on atsdr.spl's Rank column, which the DoD list has no
# equivalent of.
import duckdb
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DB_PATH = PROJECT_ROOT / 'data' / 'toxexposure.duckdb'
CHEMBL_SQLITE_PATH = PROJECT_ROOT / 'data' / 'chembl' / 'raw' / 'chembl_37.db'

SCHEMAS = ["ctd", "bindingdb", "atsdr", "dod"]


def connect(read_only: bool = False) -> duckdb.DuckDBPyConnection:
    """Open the unified project database with the ChEMBL SQLite dump attached read-only."""
    con = duckdb.connect(str(DB_PATH), read_only=read_only)
    con.execute("SET enable_progress_bar=false;")

    if not read_only:
        for schema in SCHEMAS:
            con.execute(f"CREATE SCHEMA IF NOT EXISTS {schema};")

    con.execute("INSTALL sqlite; LOAD sqlite;")
    # ATTACH is metadata-only (not a write to DB_PATH itself), so this is fine even read-only -
    # chembl is opened READ_ONLY regardless of the outer connection's mode.
    if "chembl" not in {row[0] for row in con.execute("SELECT database_name FROM duckdb_databases()").fetchall()}:
        con.execute(f"ATTACH '{CHEMBL_SQLITE_PATH.as_posix()}' AS chembl (TYPE SQLITE, READ_ONLY);")

    return con
