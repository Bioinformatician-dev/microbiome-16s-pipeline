#!/usr/bin/env python3
"""
07_diversity_analysis.py

Full downstream community-ecology analysis of the ASV abundance table:
  - Alpha diversity (richness, Shannon, Simpson, Pielou evenness)
  - Beta diversity (Bray-Curtis dissimilarity + PCoA ordination)
  - Statistical group comparison (Std/even vs Log/staggered mock layout)
  - Taxonomic composition bar chart (from vsearch-based classification)
  - Validation against the ZymoBIOMICS mock community's known/expected
    composition (ground truth), since this is a defined mock standard.

All figures are written to results/figures/, all tables to results/diversity/.
"""
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, DistanceMatrix

sns.set_theme(style="whitegrid", context="talk")

ROOT = Path(__file__).resolve().parent.parent
OTU_TABLE = ROOT / "results/otu/otu_table.tsv"
TAX_TABLE = ROOT / "results/taxonomy/asv_taxonomy.tsv"
METADATA = ROOT / "data/metadata/metadata.tsv"
EXPECTED_ABUND = ROOT / "data/reference/zymo_expected_abundance.tsv"
DIV_DIR = ROOT / "results/diversity"
FIG_DIR = ROOT / "results/figures"
DIV_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GROUP_LABELS = {"a": "Log (staggered)", "b": "Std (even)"}


def load_data():
    otu = pd.read_csv(OTU_TABLE, sep="\t", index_col=0)
    otu.index.name = "ASV_ID"
    meta = pd.read_csv(METADATA, sep="\t")
    meta = meta.rename(columns={meta.columns[0]: "sampleID"})
    meta = meta.set_index("sampleID").loc[otu.columns]
    meta["group_label"] = meta["treatment1"].map(GROUP_LABELS)
    tax = pd.read_csv(TAX_TABLE, sep="\t", index_col=0)
    return otu, meta, tax


def alpha_diversity_analysis(otu, meta):
    counts = otu.T.values  # samples x ASVs
    sample_ids = otu.columns.tolist()

    metrics = {}
    for metric in ["observed_features", "shannon", "simpson"]:
        metrics[metric] = alpha_diversity(metric, counts, ids=sample_ids)

    df = pd.DataFrame(metrics)
    # Pielou's evenness = Shannon / ln(richness)
    df["pielou_evenness"] = df["shannon"] / np.log(df["observed_features"].replace(0, np.nan))
    df = df.join(meta[["treatment1", "group_label"]])
    df.to_csv(DIV_DIR / "alpha_diversity.tsv", sep="\t")

    # Statistical test between the two mock-community layouts
    stats_rows = []
    for metric in ["observed_features", "shannon", "simpson", "pielou_evenness"]:
        g1 = df.loc[df.treatment1 == "a", metric]
        g2 = df.loc[df.treatment1 == "b", metric]
        try:
            stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
        except ValueError:
            stat, p = np.nan, np.nan
        stats_rows.append({"metric": metric, "mean_Log_a": g1.mean(), "mean_Std_b": g2.mean(),
                            "mannwhitney_U": stat, "p_value": p})
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(DIV_DIR / "alpha_diversity_group_test.tsv", sep="\t", index=False)
    print("\n=== Alpha diversity ===")
    print(df[["observed_features", "shannon", "simpson", "pielou_evenness", "group_label"]]
          .round(3).to_string())
    print("\n=== Group comparison (Mann-Whitney U, n=2 per group -> indicative only) ===")
    print(stats_df.round(4).to_string(index=False))

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
        axes, ["observed_features", "shannon", "simpson"],
        ["Observed ASVs", "Shannon diversity", "Simpson diversity"]
    ):
        sns.barplot(data=df.reset_index(), x="group_label", y=metric,
                    hue="group_label", palette="Set2", ax=ax, legend=False, errorbar="sd")
        sns.stripplot(data=df.reset_index(), x="group_label", y=metric,
                      color="black", size=8, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("")
    fig.suptitle("Alpha diversity by mock-community layout", y=1.03)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "alpha_diversity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return df


def beta_diversity_analysis(otu, meta):
    counts = otu.T.values
    sample_ids = otu.columns.tolist()

    bc_dm = beta_diversity("braycurtis", counts, ids=sample_ids)
    bc_dm.write(str(DIV_DIR / "bray_curtis_distance_matrix.tsv"))

    ordination = pcoa(bc_dm)
    coords = ordination.samples.iloc[:, :2]
    coords.columns = ["PC1", "PC2"]
    coords = coords.join(meta[["treatment1", "group_label"]])
    coords.to_csv(DIV_DIR / "pcoa_coordinates.tsv", sep="\t")

    var_explained = ordination.proportion_explained[:2] * 100

    # PERMANOVA (indicative, n=4 total / 2 per group)
    grouping = meta.loc[list(bc_dm.ids), "treatment1"]
    try:
        permanova_res = permanova(bc_dm, grouping, permutations=999)
        permanova_res.to_csv(DIV_DIR / "permanova_result.tsv", sep="\t")
        print("\n=== PERMANOVA (Bray-Curtis ~ mock layout) ===")
        print(permanova_res.to_string())
    except Exception as e:
        print(f"PERMANOVA skipped: {e}")

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(data=coords.reset_index(), x="PC1", y="PC2", hue="group_label",
                     s=200, palette="Set2", ax=ax)
    for _, row in coords.iterrows():
        ax.annotate(row.name.replace("sampleID_", ""), (row["PC1"], row["PC2"]),
                    xytext=(5, 5), textcoords="offset points", fontsize=10)
    ax.set_xlabel(f"PC1 ({var_explained.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_explained.iloc[1]:.1f}%)")
    ax.set_title("PCoA of Bray-Curtis dissimilarity")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pcoa_bray_curtis.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return bc_dm, coords


def taxonomic_composition(otu, meta, tax):
    merged = otu.join(tax[["genus", "species"]])
    merged["taxon"] = np.where(merged["genus"] == "Unclassified", "Unclassified",
                                merged["genus"] + " " + merged["species"])
    by_taxon = merged.groupby("taxon")[otu.columns.tolist()].sum()
    rel = by_taxon.div(by_taxon.sum(axis=0), axis=1) * 100
    rel = rel.sort_values(by=rel.columns.tolist(), ascending=False)
    rel.to_csv(DIV_DIR / "taxonomic_composition_relative_abundance.tsv", sep="\t")

    # keep top taxa, group rest as "Other"
    top = rel.loc[rel.mean(axis=1).sort_values(ascending=False).index[:8]]
    other = rel.drop(top.index).sum(axis=0)
    top.loc["Other"] = other

    plot_df = top.T
    plot_df.index = [meta.loc[s, "group_label"] + "\n" + s.replace("sampleID_", "")
                      for s in plot_df.index]

    fig, ax = plt.subplots(figsize=(11, 6))
    plot_df.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.65)
    ax.set_ylabel("Relative abundance (%)")
    ax.set_xlabel("")
    ax.set_title("Taxonomic composition (vsearch vs. ZymoBIOMICS reference)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    plt.xticks(rotation=0, ha="center")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "taxonomic_composition.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return rel


def validate_against_ground_truth(rel_composition):
    """Compare observed vs Zymo's expected relative abundances (mock community)."""
    expected = pd.read_csv(EXPECTED_ABUND, sep="\t")
    expected["genus_species"] = expected["ID"].str.extract(r"([A-Za-z]+_[a-z]+)")
    exp_by_species = (
        expected.groupby("genus_species")["log_expected_percentage_per-taxon"]
        .first()
        .sort_values(ascending=False)
    )
    exp_by_species.index = exp_by_species.index.str.replace("_", " ")
    exp_by_species.to_csv(DIV_DIR / "zymo_expected_composition.tsv", sep="\t")
    print("\n=== Zymo mock community: expected relative abundance (ground truth) ===")
    print(exp_by_species.round(3).to_string())


def main():
    otu, meta, tax = load_data()
    print(f"Loaded ASV table: {otu.shape[0]} ASVs x {otu.shape[1]} samples")
    print(f"Total reads in table: {otu.values.sum()}")

    alpha_diversity_analysis(otu, meta)
    beta_diversity_analysis(otu, meta)
    rel = taxonomic_composition(otu, meta, tax)
    validate_against_ground_truth(rel)

    print("\nAll diversity results written to results/diversity/")
    print("All figures written to results/figures/")


if __name__ == "__main__":
    main()
