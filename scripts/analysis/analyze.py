# analyze.py
# Author: Jagger Hershey
#
# Script for analyzing the CTD chemical-gene/chemical-disease interaction tables in the
# unified DuckDB database: field discovery, per-interaction-type deep dives, disease-linkage
# convergence, and chemical-gene-disease network graphs. Meta/whole-database exploration
# (schema summary, overall unique counts, study/chemical gene distributions) lives in meta.py.
#
# Every heavy filter/group-by runs in DuckDB against ctd.* directly; only the already-small,
# aggregated result of each query is pulled into pandas for seaborn/matplotlib.
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'common'))
from db import connect


def file_analysis(con):
    print("=====Chemical-Gene Interactions File Fields Analysis=====")

    # Discover edge types in Chemical-Gene Interactions: split the pipe-delimited
    # InteractionActions list, then each action's "prefix^suffix" (e.g. increases^activity).
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _chem_gene_actions AS
        SELECT ChemicalName, GeneSymbol, GeneForms, Organism,
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

    print("Most Common Interaction Types (Suffix):")
    print(con.execute("SELECT Suffix, count(*) AS n FROM _chem_gene_actions_split GROUP BY Suffix ORDER BY n DESC LIMIT 25").df())
    print("Unique Interaction Variation:")
    print(con.execute("SELECT DISTINCT Prefix FROM _chem_gene_actions_split").df())
    print("Unique Interaction Types:")
    print(con.execute("SELECT DISTINCT Suffix FROM _chem_gene_actions_split").df())
    print("Counts for Each Prefix:")
    print(con.execute("SELECT Prefix, count(*) AS n FROM _chem_gene_actions_split GROUP BY Prefix ORDER BY n DESC").df())

    print("Unique GeneForms Field Options:")
    print(con.execute("""
        SELECT DISTINCT UNNEST(string_split(GeneForms, '|')) AS GeneForm
        FROM ctd.chem_gene_ixns WHERE GeneForms IS NOT NULL
    """).df())
    print("Interaction Field Examples (irrelevant - concatenation of chemical+prefix+suffix+gene):")
    print(con.execute("SELECT Interaction FROM ctd.chem_gene_ixns LIMIT 5").df())

    type_analysis(con, "binding")
    type_analysis(con, "folding")

    print("=====Chemical-Disease Interactions File Fields Analysis=====")
    print("Types of Direct Evidence (therapeutic vs marker/mechanism):")
    print(con.execute("SELECT DISTINCT DirectEvidence FROM ctd.chemicals_diseases").df())
    print("Sample of InferenceGeneSymbol field (inferred relationships require a chem-gene-disease path):")
    print(con.execute("SELECT InferenceGeneSymbol FROM ctd.chemicals_diseases LIMIT 5").df())
    print("Sample InferenceScore (higher = more atypical similarity vs a scale-free random network):")
    print(con.execute("SELECT InferenceScore FROM ctd.chemicals_diseases LIMIT 5").df())

    most_common_chems = con.execute("""
        SELECT ChemicalName FROM ctd.chemicals_diseases
        GROUP BY ChemicalName ORDER BY count(*) DESC LIMIT 3
    """).df()['ChemicalName'].tolist()
    print(most_common_chems)

    chemical_analysis(con, "Benzo(a)pyrene")
    for chem in most_common_chems:
        chemical_analysis(con, chem)


def type_analysis(con, interaction_type):
    print(f"====={interaction_type.capitalize()} Interaction Type Analysis=====")
    base = "_chem_gene_actions_split"

    print(f"Number of {interaction_type.capitalize()} Interactions:",
          con.execute(f"SELECT count(*) FROM {base} WHERE Suffix = ?", [interaction_type]).fetchone()[0])
    print(f"Unique Number of {interaction_type.capitalize()} Gene Interactions:",
          con.execute(f"SELECT count(DISTINCT GeneSymbol) FROM {base} WHERE Suffix = ?", [interaction_type]).fetchone()[0])
    print("Counts for Each Prefix:")
    print(con.execute(f"SELECT Prefix, count(*) AS n FROM {base} WHERE Suffix = ? GROUP BY Prefix ORDER BY n DESC", [interaction_type]).df())
    print(f"Most Common Gene {interaction_type.capitalize()}:")
    print(con.execute(f"SELECT GeneSymbol, count(*) AS n FROM {base} WHERE Suffix = ? GROUP BY GeneSymbol ORDER BY n DESC LIMIT 5", [interaction_type]).df())
    print(f"Least Common Gene {interaction_type.capitalize()}:")
    print(con.execute(f"SELECT GeneSymbol, count(*) AS n FROM {base} WHERE Suffix = ? GROUP BY GeneSymbol ORDER BY n ASC LIMIT 5", [interaction_type]).df())
    print(f"Most Common Chemical {interaction_type.capitalize()} Interactions:")
    print(con.execute(f"SELECT ChemicalName, count(*) AS n FROM {base} WHERE Suffix = ? GROUP BY ChemicalName ORDER BY n DESC LIMIT 10", [interaction_type]).df())

    print(f"----- Filtered for Homo sapiens -----")
    human_filter = f"{base} WHERE Suffix = ? AND Organism ILIKE '%Homo sapiens%'"
    print(f"Number of Human {interaction_type.capitalize()} Interactions:",
          con.execute(f"SELECT count(*) FROM {human_filter}", [interaction_type]).fetchone()[0])
    print(f"Unique Number of Human {interaction_type.capitalize()} Gene Interactions:",
          con.execute(f"SELECT count(DISTINCT GeneSymbol) FROM {human_filter}", [interaction_type]).fetchone()[0])
    print(f"Counts for Each Prefix (Human):")
    print(con.execute(f"SELECT Prefix, count(*) AS n FROM {human_filter} GROUP BY Prefix ORDER BY n DESC", [interaction_type]).df())
    print(f"Most Common Human Gene {interaction_type.capitalize()}:")
    print(con.execute(f"SELECT GeneSymbol, count(*) AS n FROM {human_filter} GROUP BY GeneSymbol ORDER BY n DESC LIMIT 5", [interaction_type]).df())
    print(f"Least Common Human Gene {interaction_type.capitalize()}:")
    print(con.execute(f"SELECT GeneSymbol, count(*) AS n FROM {human_filter} GROUP BY GeneSymbol ORDER BY n ASC LIMIT 5", [interaction_type]).df())
    print(f"Most Common Human Chemical {interaction_type.capitalize()} Interactions:")
    print(con.execute(f"SELECT ChemicalName, count(*) AS n FROM {human_filter} GROUP BY ChemicalName ORDER BY n DESC LIMIT 10", [interaction_type]).df())


def geneform_type_distribution(con, human_only=False, top_n_forms=6):
    """
    Breaks down each CTD interaction type (Suffix) by the molecular level it was recorded at
    (GeneForms: mRNA, protein, gene, promoter, ...) - e.g. does "expression" show up almost
    entirely as mRNA-level records while "binding"/"folding" show up at the protein level?

    GeneForms is its own pipe-delimited list per row, independent of InteractionActions, so
    each row's forms are cross-joined against that same row's interaction types (a row with 2
    gene forms and 3 interaction actions contributes 6 Suffix-GeneForm pairs) - the source data
    doesn't pair individual actions to individual forms any more finely than "these forms and
    these actions both describe this record."

    Prints the overall GeneForm frequency, then a per-Suffix "signature" table (its dominant
    GeneForm and what share of that type's rows it accounts for), and draws a 100%-stacked bar
    chart per Suffix so types can be compared by composition regardless of total volume. Only
    the top_n_forms most common gene forms get their own stack segment/color; the rest are
    grouped into "other" to keep the legend readable.
    """
    organism_filter = "AND Organism ILIKE '%Homo sapiens%'" if human_only else ""
    label = "Human-Only" if human_only else "All Organisms"

    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _geneform_rows AS
        SELECT row_number() OVER () AS rid, GeneForms, InteractionActions
        FROM ctd.chem_gene_ixns
        WHERE InteractionActions IS NOT NULL AND GeneForms IS NOT NULL {organism_filter}
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _geneform_row_suffix AS
        SELECT rid, lower(trim(split_part(UNNEST(string_split(InteractionActions, '|')), '^', 2))) AS Suffix
        FROM _geneform_rows
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _geneform_row_form AS
        SELECT rid, UNNEST(string_split(GeneForms, '|')) AS GeneForm
        FROM _geneform_rows
    """)

    print(f"===== GeneForm Frequency Overall ({label}) =====")
    overall = con.execute("""
        SELECT GeneForm, count(*) AS n FROM _geneform_row_form GROUP BY GeneForm ORDER BY n DESC
    """).df()
    print(overall)

    crosstab = con.execute("""
        SELECT rs.Suffix, rf.GeneForm, count(*) AS n
        FROM _geneform_row_suffix rs JOIN _geneform_row_form rf ON rf.rid = rs.rid
        GROUP BY rs.Suffix, rf.GeneForm
    """).df()

    totals = crosstab.groupby('Suffix')['n'].sum().rename('total')
    crosstab = crosstab.merge(totals, on='Suffix')
    crosstab['pct'] = 100 * crosstab['n'] / crosstab['total']

    print(f"\n===== Interaction Type x GeneForm Signature ({label}) =====")
    print("(dominant GeneForm per Suffix, and what share of that Suffix's rows it accounts for)")
    signature = (
        crosstab.sort_values(['Suffix', 'n'], ascending=[True, False])
                .groupby('Suffix', as_index=True).first()[['total', 'GeneForm', 'pct']]
                .rename(columns={'GeneForm': 'DominantGeneForm', 'pct': 'DominantSharePct'})
                .sort_values('total', ascending=False)
    )
    print(signature)

    # Bucket all but the top_n_forms most common forms into "other" so the chart legend stays readable
    top_forms = overall['GeneForm'].head(top_n_forms).tolist()
    crosstab['GeneFormBucket'] = crosstab['GeneForm'].where(crosstab['GeneForm'].isin(top_forms), 'other')
    pivot = crosstab.groupby(['Suffix', 'GeneFormBucket'])['n'].sum().unstack(fill_value=0)
    pivot = pivot.loc[totals.sort_values(ascending=False).index]
    pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100

    ordered_cols = [f for f in top_forms if f in pivot_pct.columns] + (['other'] if 'other' in pivot_pct.columns else [])
    pivot_pct = pivot_pct[ordered_cols]

    pivot_pct.plot(kind='bar', stacked=True, figsize=(16, 7), colormap='tab20')
    plt.ylabel('Share of Interactions (%)')
    plt.xlabel('Interaction Type')
    plt.xticks(rotation=90, fontsize=8)
    plt.title(f'CTD Interaction Type Composition by GeneForm ({label})')
    plt.legend(title='GeneForm', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.show()

    return signature


def _build_gene_hub_stats(con, human_only):
    """
    Populate _all_pairs (distinct ChemicalName/GeneSymbol pairs across ALL interaction types
    combined) and _gene_hub_stats (GeneSymbol, n_chemicals - the "hub score": how many distinct
    chemicals interact with that gene anywhere in CTD).
    """
    where_clause = "WHERE Organism ILIKE '%Homo sapiens%'" if human_only else ""
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _all_pairs AS
        SELECT DISTINCT ChemicalName, GeneSymbol
        FROM ctd.chem_gene_ixns
        {where_clause}
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE _gene_hub_stats AS
        SELECT GeneSymbol, count(DISTINCT ChemicalName) AS n_chemicals
        FROM _all_pairs GROUP BY GeneSymbol
    """)


def gene_hub_distribution(con, human_only=True, hub_min_chemicals=100, specific_max_chemicals=2):
    """
    "Hub score" for a gene = how many distinct chemicals interact with it anywhere in CTD
    (across all interaction types combined). Genes with a hub score >= hub_min_chemicals are
    "hub" genes (well-trodden ground for toxicology research); genes with a hub score
    <= specific_max_chemicals are "specific" genes (only ever tied to one or two chemicals in
    the whole database - same "sparse curation vs sparse biology" caveat as the rare
    interaction-type work applies here). Defaults (100 / 2) were picked off the actual
    distribution: ~100 chemicals marks roughly the top 3% of genes, while <=2 chemicals covers
    roughly the bottom 18%.

    Prints the full hub-score distribution (so the long tail speaks for itself, not just the
    fixed cutoffs) plus the fixed-threshold hub/specific counts and top hub genes, and draws a
    rank-frequency (log-log) plot with the two thresholds marked.
    """
    _build_gene_hub_stats(con, human_only)
    label = "Human-Only" if human_only else "All Organisms"

    stats = con.execute("""
        SELECT GeneSymbol, n_chemicals FROM _gene_hub_stats ORDER BY n_chemicals DESC
    """).df()

    print(f"===== Gene Hub-Score Distribution ({label}) =====")
    print(f"Total genes: {len(stats)}")
    print(stats['n_chemicals'].describe(percentiles=[.5, .75, .9, .95, .99]))

    n_hub = (stats['n_chemicals'] >= hub_min_chemicals).sum()
    n_specific = (stats['n_chemicals'] <= specific_max_chemicals).sum()
    print(f"\nHub genes (>= {hub_min_chemicals} chemicals): {n_hub} ({n_hub / len(stats) * 100:.1f}%)")
    print(f"Specific genes (<= {specific_max_chemicals} chemicals): {n_specific} ({n_specific / len(stats) * 100:.1f}%)")

    print(f"\nTop 15 hub genes:")
    print(stats.head(15))
    print(f"\nSample of specific genes (first 15 alphabetically):")
    print(stats[stats['n_chemicals'] <= specific_max_chemicals].sort_values('GeneSymbol').head(15))

    plt.figure(figsize=(9, 6))
    plt.plot(range(1, len(stats) + 1), stats['n_chemicals'], marker='.', linestyle='none', markersize=3)
    plt.xscale('log')
    plt.yscale('log')
    plt.axhline(hub_min_chemicals, color='#2ca02c', linestyle='--', linewidth=1,
                label=f'Hub threshold (>= {hub_min_chemicals})')
    plt.axhline(specific_max_chemicals, color='#d62728', linestyle='--', linewidth=1,
                label=f'Specific threshold (<= {specific_max_chemicals})')
    plt.xlabel('Gene rank (by number of distinct chemicals interacting with it)')
    plt.ylabel('Number of distinct chemicals (log scale)')
    plt.title(f'CTD Gene Hub-Score Distribution ({label})')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return stats


def chemical_hub_specific_overlap(con, human_only=True, hub_min_chemicals=100, specific_max_chemicals=2, top_n=15):
    """
    Rather than blending a chemical's gene targets into one averaged "hub-score" (hard to read
    and not very actionable), this counts, separately, how many hub genes vs how many specific
    genes each chemical touches (see gene_hub_distribution for what hub/specific mean), and
    ranks chemicals independently on each count - "which chemicals are common among the hub
    genes" and "which chemicals are common among the specific/narrowly-studied genes" as two
    plain lists instead of one combined metric.
    """
    _build_gene_hub_stats(con, human_only)
    label = "Human-Only" if human_only else "All Organisms"

    hub_chems = con.execute("""
        SELECT ap.ChemicalName, count(DISTINCT ap.GeneSymbol) AS n_hub_genes
        FROM _all_pairs ap JOIN _gene_hub_stats ghs ON ghs.GeneSymbol = ap.GeneSymbol
        WHERE ghs.n_chemicals >= ?
        GROUP BY ap.ChemicalName ORDER BY n_hub_genes DESC LIMIT ?
    """, [hub_min_chemicals, top_n]).df()

    specific_chems = con.execute("""
        SELECT ap.ChemicalName, count(DISTINCT ap.GeneSymbol) AS n_specific_genes
        FROM _all_pairs ap JOIN _gene_hub_stats ghs ON ghs.GeneSymbol = ap.GeneSymbol
        WHERE ghs.n_chemicals <= ?
        GROUP BY ap.ChemicalName ORDER BY n_specific_genes DESC LIMIT ?
    """, [specific_max_chemicals, top_n]).df()

    print(f"===== Chemicals Most Common Among Hub Genes (>= {hub_min_chemicals} chemicals/gene) ({label}) =====")
    print(hub_chems)
    print(f"\n===== Chemicals Most Common Among Specific Genes (<= {specific_max_chemicals} chemicals/gene) ({label}) =====")
    print(specific_chems)

    _, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].barh(hub_chems['ChemicalName'][::-1], hub_chems['n_hub_genes'][::-1], color='#2ca02c')
    axes[0].set_xlabel('Number of distinct hub genes touched')
    axes[0].set_title(f'Top {top_n} Chemicals by Hub-Gene Reach')

    axes[1].barh(specific_chems['ChemicalName'][::-1], specific_chems['n_specific_genes'][::-1], color='#d62728')
    axes[1].set_xlabel('Number of distinct specific genes touched')
    axes[1].set_title(f'Top {top_n} Chemicals by Specific-Gene Reach')

    plt.suptitle(f'Chemicals Common Among Hub vs Specific Genes ({label})')
    plt.tight_layout()
    plt.show()

    return {"hub": hub_chems, "specific": specific_chems}


def _build_type_pairs(con, interaction_type, human_only):
    """Populate temp table _type_pairs(ChemicalName, GeneSymbol, Prefix) for one CTD interaction type."""
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

    organism_filter = "AND Organism ILIKE '%Homo sapiens%'" if human_only else ""
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _type_pairs AS
        SELECT DISTINCT ChemicalName, GeneSymbol, Prefix
        FROM _chem_gene_actions_split
        WHERE Suffix = ? {organism_filter}
    """, [interaction_type])


def type_disease_links(con, interaction_type, human_only=True, limit=None):
    """
    Show the concrete chemical-gene pairs for a CTD interaction type, join to
    ctd.genes_diseases (DirectEvidence only) to find diseases the gene is tied to, then
    classify each resulting (chemical, gene, disease) candidate into exactly one of three
    tiers - because CTD's own chemicals_diseases table already contains ~9.66M rows that CTD
    itself infers using this SAME chemical->gene->disease logic (DirectEvidence NULL,
    InferenceGeneSymbol/InferenceScore populated), separate from ~109K directly-curated rows.
    Reconstructing one of CTD's own inferred rows isn't a new finding, so it gets its own tier
    instead of being lumped in with genuinely new candidates:

      - direct_evidence:  the chemical ALSO has its own DIRECT evidence (marker/mechanism or
                          therapeutic) for that same disease - two independent, literature-
                          curated evidence lines agreeing. The strongest tier.
      - already_inferred: CTD's chemicals_diseases already has an INFERRED row (DirectEvidence
                          NULL) for this exact chemical/disease pair with InferenceGeneSymbol
                          equal to OUR gene - CTD's own inference engine already computed this
                          exact chemical->gene->disease connection. Reported for transparency,
                          but not a new finding - this analysis is just reconstructing something
                          already published in CTD.
      - novel:            CTD's chemicals_diseases has NEITHER a direct-evidence row NOR an
                          inferred-via-this-gene row for that chemical/disease pair - CTD's own
                          inference engine did not surface this connection either, so this is a
                          genuinely new candidate hypothesis unique to this analysis.

    NOTE: in practice "novel" is usually empty, and that's expected, not a bug. CTD's own
    inference engine fires on ANY chemical-gene relationship anywhere in chem_gene_ixns,
    regardless of interaction type - so if a chemical and gene co-occur through even one of the
    other 52 interaction types (very common for well-studied chemicals), CTD has already
    inferred every disease that gene has direct evidence for. Since this analysis's candidates
    are always a same-or-smaller subset of that full relationship, CTD's inference has
    essentially always gotten there first. The real value here is the "direct_evidence" vs
    "already_inferred" split (how much of what we'd otherwise call a "finding" is actually
    brand new evidence-agreement vs. CTD having already computed the same thing) - not
    discovering links CTD doesn't already know about.

    Pass limit=N to cap the printed/returned tables for high-volume types (e.g. binding);
    the summary counts always reflect the true unfiltered totals.
    """
    _build_type_pairs(con, interaction_type, human_only)

    label = "Human-Only" if human_only else "All Organisms"
    limit_clause = f"LIMIT {int(limit)}" if limit else ""

    print(f"====={interaction_type.capitalize()} Chemical-Gene Pairs ({label})=====")
    n_pairs = con.execute("SELECT count(*) FROM _type_pairs").fetchone()[0]
    print(f"Unique chemical-gene pairs: {n_pairs}")
    print(con.execute(f"SELECT * FROM _type_pairs ORDER BY ChemicalName, GeneSymbol {limit_clause}").df())

    print(f"\n====={interaction_type.capitalize()} Gene Disease Associations (Direct Evidence Only)=====")
    n_gene_disease = con.execute("""
        SELECT count(*) FROM _type_pairs tp
        JOIN ctd.genes_diseases gd ON gd.GeneSymbol = tp.GeneSymbol
        WHERE gd.DirectEvidence IS NOT NULL
    """).fetchone()[0]
    print(f"Gene-disease direct-evidence rows found: {n_gene_disease}")
    gene_disease = con.execute(f"""
        SELECT tp.ChemicalName, tp.GeneSymbol, tp.Prefix, gd.DiseaseName, gd.DirectEvidence AS GeneEvidence
        FROM _type_pairs tp
        JOIN ctd.genes_diseases gd ON gd.GeneSymbol = tp.GeneSymbol
        WHERE gd.DirectEvidence IS NOT NULL
        ORDER BY tp.ChemicalName, tp.GeneSymbol
        {limit_clause}
    """).df()
    print(gene_disease)

    classified = con.execute("""
        SELECT tp.ChemicalName, tp.GeneSymbol, tp.Prefix, gd.DiseaseName,
               gd.DirectEvidence AS GeneEvidence,
               cd_direct.DirectEvidence AS ChemDirectEvidence,
               cd_inferred.InferenceScore AS CtdInferenceScore,
               CASE
                   WHEN cd_direct.DiseaseName IS NOT NULL THEN 'direct_evidence'
                   WHEN cd_inferred.DiseaseName IS NOT NULL THEN 'already_inferred'
                   ELSE 'novel'
               END AS Category
        FROM _type_pairs tp
        JOIN ctd.genes_diseases gd ON gd.GeneSymbol = tp.GeneSymbol AND gd.DirectEvidence IS NOT NULL
        LEFT JOIN ctd.chemicals_diseases cd_direct
               ON cd_direct.ChemicalName = tp.ChemicalName AND cd_direct.DiseaseName = gd.DiseaseName
              AND cd_direct.DirectEvidence IS NOT NULL
        LEFT JOIN ctd.chemicals_diseases cd_inferred
               ON cd_inferred.ChemicalName = tp.ChemicalName AND cd_inferred.DiseaseName = gd.DiseaseName
              AND cd_inferred.InferenceGeneSymbol = tp.GeneSymbol
    """).df()

    counts = classified['Category'].value_counts()
    print(f"\n=====Chemical<->Gene<->Disease Candidates for {interaction_type.capitalize()}, Classified=====")
    print(f"Direct-evidence corroborated (strongest): {counts.get('direct_evidence', 0)}")
    print(f"Already CTD-inferred (not novel, reconstructs an existing CTD row): {counts.get('already_inferred', 0)}")
    print(f"Novel (unique to this analysis - not in CTD in any form): {counts.get('novel', 0)}")

    direct_evidence = classified[classified['Category'] == 'direct_evidence'].sort_values(['ChemicalName', 'GeneSymbol'])
    already_inferred = classified[classified['Category'] == 'already_inferred'].sort_values(['ChemicalName', 'GeneSymbol'])
    novel = classified[classified['Category'] == 'novel'].sort_values(['ChemicalName', 'GeneSymbol'])

    print(f"\n--- Direct-evidence corroborated ---")
    print(direct_evidence.head(limit) if limit else direct_evidence)

    print(f"\n--- Already CTD-inferred (filtered out of \"novel\" below - shown for transparency) ---")
    print(already_inferred.head(limit) if limit else already_inferred)

    print(f"\n--- Novel: unique to this analysis ---")
    print(novel.head(limit) if limit else novel)

    return {"pairs": n_pairs, "gene_disease": gene_disease,
            "direct_evidence": direct_evidence, "already_inferred": already_inferred, "novel": novel}


def build_type_disease_network(con, interaction_type, human_only=True, limit=25, diseases_per_gene=3):
    """
    Chemical -> gene -> disease network graph for one CTD interaction type.

    Edges:
      chemical -> gene    the interaction itself, colored by direction (increases/decreases/affects)
      gene -> disease     direct evidence only (marker/mechanism or therapeutic) - always drawn
      chemical -> disease drawn in one of two ways (see type_disease_links for the full
                          three-tier explanation of these categories):
                            - orange solid: "direct_evidence" - the chemical ALSO has its own
                              independent direct evidence for that disease. Strongest tier.
                            - purple dotted: "novel" - CTD has NEITHER a direct-evidence row NOR
                              an inferred-via-this-gene row for that chemical/disease pair, i.e.
                              this triangle isn't already published anywhere in CTD.
                          "already_inferred" candidates (CTD's own inference engine already
                          computed this exact chemical->gene->disease link) get NO extra edge at
                          all - filtered out rather than drawn, since presenting them as a
                          "finding" would just be reconstructing something already in CTD.

    RANKING (how `limit` decides which chemical-gene pairs make the cut):
    High-volume types (binding, folding, ...) can have tens of thousands of chemical-gene pairs -
    far too many to plot - so every pair is scored and only the top `limit` are drawn. Each pair
    is scored on two numbers, computed by joining it out to every disease its gene has direct
    evidence for and classifying each resulting triad (same three tiers as above):
      1. n_direct_evidence - how many diseases this pair corroborates via independent direct
         evidence. This is the primary sort key: pairs with real, independently-curated
         corroboration are the most trustworthy and get drawn first.
      2. n_novel - how many diseases this pair surfaces that AREN'T already in CTD in any form.
         This is the tiebreaker: among pairs with equal (often zero) direct-evidence support,
         pairs that surface more genuinely new candidate disease links are more interesting to
         look at than pairs that just repeat what CTD's own inference engine already knows.
    n_already_inferred deliberately does NOT contribute to the ranking score at all - a pair that
    only reconstructs CTD's existing inferred rows shouldn't outrank a pair with real direct
    evidence or real novelty, no matter how many diseases it touches.

    `diseases_per_gene` caps how many disease nodes a single promiscuous gene can add to the
    drawing - direct_evidence diseases are kept first, novel diseases next, and already_inferred
    diseases fill any remaining slots last (so if a gene has to be trimmed, the already-inferred
    ones are dropped before anything more meaningful), keeping a few "hub" genes from swamping
    the layout.
    """
    _build_type_pairs(con, interaction_type, human_only)
    label = "Human-Only" if human_only else "All Organisms"

    ranked_pairs = con.execute("""
        SELECT tp.ChemicalName, tp.GeneSymbol, tp.Prefix,
               count(DISTINCT CASE WHEN cd_direct.DiseaseName IS NOT NULL THEN gd.DiseaseName END) AS n_direct_evidence,
               count(DISTINCT CASE WHEN cd_direct.DiseaseName IS NULL AND cd_inferred.DiseaseName IS NULL
                                    THEN gd.DiseaseName END) AS n_novel
        FROM _type_pairs tp
        LEFT JOIN ctd.genes_diseases gd ON gd.GeneSymbol = tp.GeneSymbol AND gd.DirectEvidence IS NOT NULL
        LEFT JOIN ctd.chemicals_diseases cd_direct
               ON cd_direct.ChemicalName = tp.ChemicalName AND cd_direct.DiseaseName = gd.DiseaseName
              AND cd_direct.DirectEvidence IS NOT NULL
        LEFT JOIN ctd.chemicals_diseases cd_inferred
               ON cd_inferred.ChemicalName = tp.ChemicalName AND cd_inferred.DiseaseName = gd.DiseaseName
              AND cd_inferred.InferenceGeneSymbol = tp.GeneSymbol
        GROUP BY tp.ChemicalName, tp.GeneSymbol, tp.Prefix
        ORDER BY n_direct_evidence DESC, n_novel DESC, tp.ChemicalName, tp.GeneSymbol
    """).df()

    total_pairs = len(ranked_pairs)
    if limit:
        ranked_pairs = ranked_pairs.head(limit).reset_index(drop=True)

    print(f"===== Building {interaction_type.capitalize()} Network ({label}) =====")
    print(f"Chemical-gene pairs drawn: {len(ranked_pairs)} (of {total_pairs} total)")
    print("Ranked by: # direct-evidence-corroborated diseases first, then # novel diseases as "
          "tiebreaker (already-CTD-inferred diseases don't count toward the ranking at all)")

    con.register('_network_pairs_df', ranked_pairs[['ChemicalName', 'GeneSymbol', 'Prefix']])
    gene_disease = con.execute("""
        SELECT np.ChemicalName, np.GeneSymbol, gd.DiseaseName,
               CASE
                   WHEN cd_direct.DiseaseName IS NOT NULL THEN 'direct_evidence'
                   WHEN cd_inferred.DiseaseName IS NOT NULL THEN 'already_inferred'
                   ELSE 'novel'
               END AS Category
        FROM _network_pairs_df np
        JOIN ctd.genes_diseases gd ON gd.GeneSymbol = np.GeneSymbol AND gd.DirectEvidence IS NOT NULL
        LEFT JOIN ctd.chemicals_diseases cd_direct
               ON cd_direct.ChemicalName = np.ChemicalName AND cd_direct.DiseaseName = gd.DiseaseName
              AND cd_direct.DirectEvidence IS NOT NULL
        LEFT JOIN ctd.chemicals_diseases cd_inferred
               ON cd_inferred.ChemicalName = np.ChemicalName AND cd_inferred.DiseaseName = gd.DiseaseName
              AND cd_inferred.InferenceGeneSymbol = np.GeneSymbol
        ORDER BY np.GeneSymbol,
                 CASE Category WHEN 'direct_evidence' THEN 0 WHEN 'novel' THEN 1 ELSE 2 END
    """).df()
    con.unregister('_network_pairs_df')

    # Cap disease nodes per gene, keeping direct_evidence first, novel second, already_inferred
    # last (so already_inferred is what gets dropped if a hub gene has to be trimmed).
    gene_disease = gene_disease.groupby('GeneSymbol', group_keys=False).head(diseases_per_gene)
    cat_counts = gene_disease['Category'].value_counts()
    print(f"Disease nodes drawn: {gene_disease['DiseaseName'].nunique()} "
          f"({cat_counts.get('direct_evidence', 0)} direct-evidence edges, "
          f"{cat_counts.get('novel', 0)} novel-hypothesis edges, "
          f"{cat_counts.get('already_inferred', 0)} already-inferred [no chemical-disease edge drawn])")

    G = nx.Graph()
    node_kind = {}
    prefix_color = {"increases": "#2ca02c", "decreases": "#d62728", "affects": "#7f7f7f"}

    for _, row in ranked_pairs.iterrows():
        node_kind[row['ChemicalName']] = 'chemical'
        node_kind[row['GeneSymbol']] = 'gene'
        G.add_edge(row['ChemicalName'], row['GeneSymbol'],
                   color=prefix_color.get(row['Prefix'], "#7f7f7f"), style='solid', width=1.5)

    chems_by_gene = ranked_pairs.groupby('GeneSymbol')['ChemicalName'].apply(list)
    for _, row in gene_disease.iterrows():
        node_kind[row['DiseaseName']] = 'disease'
        G.add_edge(row['GeneSymbol'], row['DiseaseName'], color="#1f77b4", style='dashed', width=1.0)
        if row['Category'] == 'direct_evidence':
            for chem in chems_by_gene.get(row['GeneSymbol'], []):
                G.add_edge(chem, row['DiseaseName'], color="#ff7f0e", style='solid', width=2.5)
        elif row['Category'] == 'novel':
            for chem in chems_by_gene.get(row['GeneSymbol'], []):
                G.add_edge(chem, row['DiseaseName'], color="#9467bd", style='dotted', width=2.0)
        # already_inferred: no chemical-disease edge - filtered out, not a new finding

    kind_color = {'chemical': '#4C72B0', 'gene': '#55A868', 'disease': '#C44E52'}
    pos = nx.spring_layout(G, seed=42, k=0.6)
    plt.figure(figsize=(16, 12))

    for kind, color in kind_color.items():
        nodes = [n for n in G.nodes() if node_kind[n] == kind]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=color, node_size=400, label=kind.capitalize())

    for style in ('solid', 'dashed', 'dotted'):
        edges = [(u, v) for u, v in G.edges() if G[u][v]['style'] == style]
        if edges:
            nx.draw_networkx_edges(G, pos, edgelist=edges, style=style,
                                    edge_color=[G[u][v]['color'] for u, v in edges],
                                    width=[G[u][v]['width'] for u, v in edges])

    nx.draw_networkx_labels(G, pos, font_size=7)
    plt.title(f"{interaction_type.capitalize()} network: chemical -> gene -> disease ({label})\n"
              f"orange = direct-evidence corroborated, purple dotted = novel (not yet in CTD)")

    node_handles, _ = plt.gca().get_legend_handles_labels()
    edge_handles = [
        Line2D([0], [0], color=prefix_color['increases'], lw=1.5, label='Increases'),
        Line2D([0], [0], color=prefix_color['decreases'], lw=1.5, label='Decreases'),
        Line2D([0], [0], color=prefix_color['affects'], lw=1.5, label='Affects'),
        Line2D([0], [0], color="#1f77b4", lw=1.0, linestyle='dashed', label='Gene-disease link'),
        Line2D([0], [0], color="#ff7f0e", lw=2.5, label='Direct-evidence corroboration'),
        Line2D([0], [0], color="#9467bd", lw=2.0, linestyle='dotted', label='Novel (not in CTD)'),
    ]
    plt.legend(handles=node_handles + edge_handles, scatterpoints=1, loc='best', fontsize=8)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

    return G


def chemical_analysis(con, chemical_name):
    print(f"====={chemical_name}=====")
    print(f"Number of {chemical_name} Interactions:",
          con.execute("SELECT count(*) FROM ctd.chemicals_diseases WHERE ChemicalName = ?", [chemical_name]).fetchone()[0])
    print(f"Unique Number of {chemical_name} Disease Interactions:",
          con.execute("SELECT count(DISTINCT DiseaseName) FROM ctd.chemicals_diseases WHERE ChemicalName = ?", [chemical_name]).fetchone()[0])
    print("Counts for Each Type of Direct Evidence:")
    print(con.execute("SELECT DirectEvidence, count(*) AS n FROM ctd.chemicals_diseases WHERE ChemicalName = ? GROUP BY DirectEvidence ORDER BY n DESC", [chemical_name]).df())
    print(f"Most Common {chemical_name} Associated Disease:")
    print(con.execute("SELECT DiseaseName, count(*) AS n FROM ctd.chemicals_diseases WHERE ChemicalName = ? GROUP BY DiseaseName ORDER BY n DESC LIMIT 5", [chemical_name]).df())
    print(f"Least Common {chemical_name} Associated Disease:")
    print(con.execute("SELECT DiseaseName, count(*) AS n FROM ctd.chemicals_diseases WHERE ChemicalName = ? GROUP BY DiseaseName ORDER BY n ASC LIMIT 5", [chemical_name]).df())


def atsdr_top_chemical_patterns(con, top_n=15, human_only=True, min_shared=2, display_n=15):
    """
    Take the ATSDR Substance Priority List's top_n most toxic substances (official current-year
    Rank), match them case-insensitively to ctd.chem_gene_ixns.ChemicalName, and look for
    convergent biology across them:
      - genes hit by >= min_shared of these top substances, broken down by which interaction
        type (Suffix) each contributing substance hits that gene through
      - diseases >= min_shared of these substances are independently linked to
        (ctd.chemicals_diseases, DirectEvidence only)

    A substance with no exact ChemicalName match in CTD is reported as a coverage gap rather
    than silently dropped or fuzzy-matched to a different (possibly wrong) chemical - e.g. CTD
    may only have a specific positional isomer under a different exact name.
    """
    label = "Human-Only" if human_only else "All Organisms"
    organism_filter = "AND Organism ILIKE '%Homo sapiens%'" if human_only else ""

    atsdr_list = con.execute("""
        SELECT Rank, "Substance Name" AS Substance FROM atsdr.spl WHERE Rank <= ? ORDER BY Rank
    """, [top_n]).df()
    con.register('_atsdr_top', atsdr_list)

    con.execute("""
        CREATE OR REPLACE TEMP TABLE _atsdr_matched AS
        SELECT DISTINCT a.Rank, a.Substance, cg.ChemicalName
        FROM _atsdr_top a JOIN ctd.chem_gene_ixns cg ON lower(cg.ChemicalName) = lower(a.Substance)
    """)
    con.unregister('_atsdr_top')

    matched = con.execute("SELECT Rank, Substance FROM _atsdr_matched ORDER BY Rank").df()
    unmatched = atsdr_list[~atsdr_list['Substance'].isin(matched['Substance'])]

    print(f"===== ATSDR Top {top_n} Substances -> CTD Chemical Match ({label}) =====")
    print(f"Matched {matched['Substance'].nunique()}/{top_n} to a CTD chemical name")
    if len(unmatched):
        print("Unmatched (no exact ChemicalName match found in CTD - coverage gap, not analyzed further):")
        print(unmatched)

    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _atsdr_chem_gene_split AS
        SELECT am.Rank, am.Substance, cg.GeneSymbol,
               lower(trim(split_part(UNNEST(string_split(cg.InteractionActions, '|')), '^', 2))) AS Suffix
        FROM ctd.chem_gene_ixns cg
        JOIN _atsdr_matched am ON am.ChemicalName = cg.ChemicalName
        WHERE cg.InteractionActions IS NOT NULL {organism_filter}
    """)
    gene_share = con.execute("""
        SELECT Substance, GeneSymbol, Suffix, count(*) AS n
        FROM _atsdr_chem_gene_split GROUP BY Substance, GeneSymbol, Suffix
    """).df()

    n_substances_per_gene = gene_share.groupby('GeneSymbol')['Substance'].nunique().rename('n_substances')
    shared_genes = n_substances_per_gene[n_substances_per_gene >= min_shared].sort_values(ascending=False)

    print(f"\n===== Genes Shared Across Top {top_n} Substances ({label}) =====")
    print(f"Genes hit by >= {min_shared} of the matched substances: {len(shared_genes)}")
    for gene, n_sub in shared_genes.head(display_n).items():
        print(f"\n{gene} (hit by {n_sub} substances):")
        gene_rows = gene_share[gene_share['GeneSymbol'] == gene]
        for substance, grp in gene_rows.groupby('Substance'):
            types = ', '.join(sorted(grp['Suffix'].unique()))
            print(f"    {substance}: {types}")

    disease_pairs = con.execute("""
        SELECT DISTINCT am.Substance, cd.DiseaseName
        FROM ctd.chemicals_diseases cd
        JOIN _atsdr_matched am ON am.ChemicalName = cd.ChemicalName
        WHERE cd.DirectEvidence IS NOT NULL
    """).df()
    n_substances_per_disease = disease_pairs.groupby('DiseaseName')['Substance'].nunique().rename('n_substances')
    shared_diseases = n_substances_per_disease[n_substances_per_disease >= min_shared].sort_values(ascending=False)

    print(f"\n===== Diseases Shared Across Top {top_n} Substances ({label}) =====")
    print(f"Diseases linked (direct evidence) to >= {min_shared} of the matched substances: {len(shared_diseases)}")
    for disease, n_sub in shared_diseases.head(display_n).items():
        subs = sorted(disease_pairs.loc[disease_pairs['DiseaseName'] == disease, 'Substance'].unique())
        print(f"  {disease} ({n_sub}): {', '.join(subs)}")

    return {"matched": matched, "unmatched": unmatched, "shared_genes": shared_genes, "shared_diseases": shared_diseases}


if __name__ == "__main__":
    con = connect()
    # file_analysis(con)
    # geneform_type_distribution(con, human_only=True)
    gene_hub_distribution(con)
    chemical_hub_specific_overlap(con)
    atsdr_top_chemical_patterns(con)
    # for rare_type in ("myristoylation", "polymerization", "ribosylation"):
    #     type_disease_links(con, rare_type)
    #     build_type_disease_network(con, rare_type)
    # type_disease_links(con, "binding", limit=25)
    # build_type_disease_network(con, "binding", limit=25)
    # type_disease_links(con, "folding", limit=25)
    # build_type_disease_network(con, "folding", limit=25)
    con.close()
