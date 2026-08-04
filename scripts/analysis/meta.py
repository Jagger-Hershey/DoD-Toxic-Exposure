# meta.py
# Author: Jagger Hershey
#
# Meta/one-off exploration of the unified DuckDB database: schema summaries, overall
# unique-entity counts, and study/chemical-level gene distributions. Interaction-type-level
# CTD analysis (chemical-gene fields, disease links, networks) lives in analyze.py.
import sys
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'common'))
from db import connect


def save_schema(con, out_path=None):
    """Write a quick per-table column/type summary for every attached database."""
    if out_path is None:
        out_path = Path(__file__).resolve().parent.parent.parent / 'data' / 'schema_summary.txt'

    schema_df = con.execute("""
        SELECT table_catalog, table_schema, table_name, column_name, data_type
        FROM information_schema.columns
        WHERE table_schema NOT IN ('information_schema', 'pg_catalog')
        ORDER BY table_catalog, table_schema, table_name, ordinal_position
    """).df()

    lines = ["=====Database Schema Summary====="]
    for (catalog, schema, table), group in schema_df.groupby(
        ['table_catalog', 'table_schema', 'table_name'], sort=False
    ):
        lines.append(f"\n{catalog}.{schema}.{table} ({len(group)} columns)")
        for _, row in group.iterrows():
            lines.append(f"  {row['column_name']}: {row['data_type']}")

    summary = "\n".join(lines)
    print(summary)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(summary + "\n", encoding="utf-8")
    print(f"\nSchema summary written to {out_path}")


def summarize_bulk(con):
    print("=====Unique Counts Summary=====")
    print("Unique diseases:", con.execute("SELECT count(DISTINCT DiseaseName) FROM ctd.chemicals_diseases").fetchone()[0])
    print("Unique genes:", con.execute("SELECT count(DISTINCT GeneSymbol) FROM ctd.genes_diseases").fetchone()[0])
    print("Unique chemicals:", con.execute("SELECT count(DISTINCT ChemicalName) FROM ctd.chemicals_diseases").fetchone()[0])

    print("\n=====Common Diseases and Associated Interactions=====")
    top_disease = con.execute("""
        SELECT DiseaseName, count(*) AS n FROM ctd.genes_diseases
        GROUP BY DiseaseName ORDER BY n DESC LIMIT 25
    """).df()
    print(top_disease)
    sns.barplot(data=top_disease, x='DiseaseName', y='n')
    plt.xticks(rotation=90, fontsize=8)
    plt.ticklabel_format(style="plain", axis='y')
    plt.show()

    # Materialize the top-25 disease/gene sets once so later queries can join against them
    # instead of recomputing the same top-N aggregation repeatedly.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _top_diseases AS
        SELECT DiseaseName FROM (
            SELECT DiseaseName, count(*) AS n FROM ctd.genes_diseases
            GROUP BY DiseaseName ORDER BY n DESC LIMIT 25
        )
    """)

    print("\n=====Genes Associated With the Most Common Diseases=====")
    top_genes = con.execute("""
        SELECT gd.GeneSymbol, count(*) AS n
        FROM ctd.genes_diseases gd JOIN _top_diseases td ON td.DiseaseName = gd.DiseaseName
        GROUP BY gd.GeneSymbol ORDER BY n DESC LIMIT 25
    """).df()
    print(top_genes)

    con.execute("""
        CREATE OR REPLACE TEMP TABLE _top_genes AS
        SELECT GeneSymbol FROM (
            SELECT gd.GeneSymbol, count(*) AS n
            FROM ctd.genes_diseases gd JOIN _top_diseases td ON td.DiseaseName = gd.DiseaseName
            GROUP BY gd.GeneSymbol ORDER BY n DESC LIMIT 25
        )
    """)

    print("\n=====Pathways for Those Genes=====")
    top_pathways = con.execute("""
        SELECT gp.PathwayName, count(*) AS n
        FROM ctd.genes_pathways gp JOIN _top_genes tg ON tg.GeneSymbol = gp.GeneSymbol
        GROUP BY gp.PathwayName ORDER BY n DESC LIMIT 25
    """).df()
    print(top_pathways)

    print("\n=====Chemicals Associated With the Most Common Diseases=====")
    top_chemicals = con.execute("""
        SELECT cd.ChemicalName, count(*) AS n
        FROM ctd.chemicals_diseases cd JOIN _top_diseases td ON td.DiseaseName = cd.DiseaseName
        GROUP BY cd.ChemicalName ORDER BY n DESC LIMIT 25
    """).df()
    print(top_chemicals)


def gene_distribution(con):
    print("=====Distribution of Genes Associated Per Study and Chemical filtered by Homo Sapiens=====")
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _human_chem_gene AS
        SELECT ChemicalName, GeneSymbol, UNNEST(string_split(PubMedIDs, '|')) AS PubMedID
        FROM ctd.chem_gene_ixns
        WHERE Organism ILIKE '%Homo sapiens%' AND PubMedIDs IS NOT NULL
    """)

    sample = 10
    num_genes_per_study = con.execute("""
        SELECT PubMedID, count(DISTINCT GeneSymbol) AS GeneCount
        FROM _human_chem_gene GROUP BY PubMedID ORDER BY GeneCount DESC LIMIT ?
    """, [sample]).df()
    print(num_genes_per_study)
    sns.barplot(data=num_genes_per_study, x='PubMedID', y='GeneCount')
    plt.title(f'Top {sample} Studies with Most Genes')
    plt.xticks(rotation=90)
    plt.ticklabel_format(style="plain", axis='y')
    plt.show()

    num_genes_per_chem = con.execute("""
        SELECT ChemicalName, count(DISTINCT GeneSymbol) AS GeneCount
        FROM _human_chem_gene GROUP BY ChemicalName ORDER BY GeneCount DESC LIMIT ?
    """, [sample]).df()
    print(num_genes_per_chem)
    sns.barplot(data=num_genes_per_chem, x='ChemicalName', y='GeneCount')
    plt.title(f'Top {sample} Chemicals with Most Genes')
    plt.xticks(rotation=90)
    plt.ticklabel_format(style="plain", axis='y')
    plt.show()

    threshold = 1
    num_under_threshold = con.execute("""
        SELECT PubMedID, count(DISTINCT GeneSymbol) AS GeneCount
        FROM _human_chem_gene GROUP BY PubMedID
        HAVING count(DISTINCT GeneSymbol) <= ?
        ORDER BY GeneCount ASC
    """, [threshold]).df()
    print(f"Number of studies that focused on less than or equal to {threshold} gene(s) filtered by Homo Sapiens: {len(num_under_threshold)}")
    print(num_under_threshold.head(sample))


def interaction_type_distribution(con, human_only=False):
    """Bar chart of total CTD chemical-gene interactions for each interaction type (Suffix)."""
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _chem_gene_actions AS
        SELECT ChemicalName, GeneSymbol, Organism,
               UNNEST(string_split(InteractionActions, '|')) AS InteractionAction
        FROM ctd.chem_gene_ixns
        WHERE InteractionActions IS NOT NULL
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _chem_gene_actions_split AS
        SELECT *,
               lower(trim(split_part(InteractionAction, '^', 1))) AS Prefix,
               lower(trim(split_part(InteractionAction, '^', 2))) AS Suffix
        FROM _chem_gene_actions
    """)

    where_clause = "WHERE Organism ILIKE '%Homo sapiens%'" if human_only else ""
    label = "Human-Only" if human_only else "All Organisms"

    counts = con.execute(f"""
        SELECT Suffix, count(*) AS n
        FROM _chem_gene_actions_split
        {where_clause}
        GROUP BY Suffix ORDER BY n DESC
    """).df()

    print(f"===== Interaction Type Distribution ({label}) =====")
    print(counts)

    plt.figure(figsize=(14, 6))
    sns.barplot(data=counts, x='Suffix', y='n', color="#4C72B0")
    plt.yscale('log')
    plt.xticks(rotation=90, fontsize=8)
    plt.ylabel('Interaction Count (log scale)')
    plt.xlabel('Interaction Type')
    plt.title(f'CTD Chemical-Gene Interaction Type Distribution ({label})')
    plt.tight_layout()
    plt.show()

    return counts


if __name__ == "__main__":
    con = connect()
    # save_schema(con)
    # summarize_bulk(con)
    gene_distribution(con)
    # interaction_type_distribution(con, human_only=True)
    con.close()
