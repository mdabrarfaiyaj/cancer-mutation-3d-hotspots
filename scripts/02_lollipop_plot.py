"""
02_lollipop_plot.py
-------------------
Make the classic cancer-genomics "lollipop" plot for TP53, KRAS and PIK3CA:
  x-axis = protein position, y-axis = number of patients with a missense
  mutation there, protein domains drawn underneath.

Also makes a small bar chart: which cancer types carry mutations in each gene.

Inputs  (from script 01):
  data/processed/<GENE>_missense_counts.csv
  data/raw/<GENE>_mutations_raw.csv
Outputs:
  results/figures/<GENE>_lollipop.png / .pdf
  results/figures/all_genes_lollipop.png
  results/figures/cancer_types_per_gene.png

Run from the project folder:
  python scripts/02_lollipop_plot.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

PROC_DIR = Path("data/processed")
RAW_DIR = Path("data/raw")
FIG_DIR = Path("results/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Protein length and domain boundaries (approximate, based on UniProt annotations).
# Cite UniProt (P04637, P01116, P42336) in your README.
PROTEINS = {
    "TP53": {
        "length": 393,
        "domains": [
            ("Transactivation", 1, 61, "#9ecae1"),
            ("Pro-rich", 64, 92, "#c7e9c0"),
            ("DNA-binding domain", 94, 292, "#fdae6b"),
            ("Tetramer.", 325, 356, "#bcbddc"),
        ],
    },
    "KRAS": {
        "length": 189,
        "domains": [
            ("P-loop", 10, 17, "#fdae6b"),
            ("Switch I", 30, 38, "#9ecae1"),
            ("Switch II", 59, 76, "#c7e9c0"),
            ("HVR", 167, 189, "#bcbddc"),
        ],
    },
    "PIK3CA": {
        "length": 1068,
        "domains": [
            ("ABD", 16, 105, "#bcbddc"),
            ("RBD", 187, 289, "#c7e9c0"),
            ("C2", 330, 487, "#d9d9d9"),
            ("Helical", 517, 694, "#9ecae1"),
            ("Kinase", 797, 1068, "#fdae6b"),
        ],
    },
}

N_LABELS = 5           # how many hotspots to label per gene
STEM_COLOR = "#bdbdbd"
HEAD_COLOR = "#c0392b"


def draw_lollipop(ax, gene, counts):
    info = PROTEINS[gene]
    length = info["length"]
    ymax = counts["n_patients"].max()

    # stems + heads
    ax.vlines(counts["position"], 0, counts["n_patients"], color=STEM_COLOR, lw=0.8, zorder=1)
    ax.scatter(counts["position"], counts["n_patients"], s=18, color=HEAD_COLOR,
               edgecolor="white", linewidth=0.4, zorder=2)

    # label top hotspots
    top = counts.nlargest(N_LABELS, "n_patients")
    for _, r in top.iterrows():
        ax.annotate(f"{r['top_change']}\n(n={r['n_patients']})",
                    (r["position"], r["n_patients"]),
                    xytext=(0, 6), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.5)

    # protein bar + domains drawn below the x-axis
    bar_y, bar_h = -0.13 * ymax, 0.07 * ymax
    ax.add_patch(Rectangle((1, bar_y), length, bar_h, color="#f0f0f0", ec="#969696", lw=0.6))
    for name, start, end, color in info["domains"]:
        ax.add_patch(Rectangle((start, bar_y - bar_h * 0.25), end - start, bar_h * 1.5,
                               color=color, ec="#636363", lw=0.6))
        ax.text((start + end) / 2, bar_y - bar_h * 0.9, name,
                ha="center", va="top", fontsize=7)

    ax.set_xlim(0, length + 1)
    ax.set_ylim(-0.30 * ymax, ymax * 1.30)
    ax.set_ylabel("Patients with\nmissense mutation")
    ax.set_title(f"{gene}  ({counts['n_patients'].sum()} missense mutations, "
                 f"{len(counts)} positions)", loc="left", fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_bounds(0, ymax)
    ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= ymax])


def lollipop_figures():
    fig_all, axes = plt.subplots(3, 1, figsize=(11, 10))
    for ax_all, gene in zip(axes, PROTEINS):
        counts = pd.read_csv(PROC_DIR / f"{gene}_missense_counts.csv")

        # single-gene figure
        fig, ax = plt.subplots(figsize=(11, 3.8))
        draw_lollipop(ax, gene, counts)
        ax.set_xlabel("Protein position (amino acid)")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"{gene}_lollipop.{ext}", dpi=300)
        plt.close(fig)

        # panel in combined figure
        draw_lollipop(ax_all, gene, counts)

    axes[-1].set_xlabel("Protein position (amino acid)")
    fig_all.suptitle("Missense mutation hotspots across 32 TCGA PanCancer Atlas cancer types",
                     fontsize=12, x=0.01, ha="left")
    fig_all.tight_layout()
    fig_all.savefig(FIG_DIR / "all_genes_lollipop.png", dpi=300)
    plt.close(fig_all)
    print("Saved lollipop plots")


def cancer_type_figure(top_n=8):
    """Which cancer types carry missense mutations in each gene (patients counted once)."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, gene in zip(axes, PROTEINS):
        raw = pd.read_csv(RAW_DIR / f"{gene}_mutations_raw.csv")
        mis = raw[raw["mutation_type"] == "Missense_Mutation"]
        per_type = (mis.drop_duplicates(["patient"])
                       .groupby("cancer_type")["patient"].nunique()
                       .nlargest(top_n).sort_values())
        ax.barh(per_type.index.str.upper(), per_type.values, color="#2c7fb8")
        ax.set_title(gene, loc="left", fontweight="bold")
        ax.set_xlabel("Patients")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Top {top_n} cancer types by patients with a missense mutation",
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cancer_types_per_gene.png", dpi=300)
    plt.close(fig)
    print("Saved cancer-type bar chart")


if __name__ == "__main__":
    lollipop_figures()
    cancer_type_figure()
    print(f"Done. Figures are in {FIG_DIR}/")
