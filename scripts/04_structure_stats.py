"""
04_structure_stats.py
---------------------
Two simple structural tests, using the per-residue tables from script 03.

Test 1 - Do hotspots cluster in 3D?
  Take the top-K most mutated residues and compute their average pairwise
  CA-CA distance. Compare with 10,000 random sets of K residues from the
  same structure (permutation test). Smaller than random = clustered.
  For the AlphaFold model (PIK3CA), only confidently predicted residues
  (pLDDT >= 70) are used.

Test 2 - Are hotspots close to the protein's functional partner?
  TP53: distance to DNA.   KRAS: distance to GDP (the nucleotide pocket).
  Each hotspot is called "contact" if any atom is within 5 A of the partner.

Outputs:
  results/clustering_test.csv
  results/distance_to_partner.csv
  results/figures/clustering_test.png
  results/figures/distance_to_partner.png

Run from the project folder:
  python scripts/04_structure_stats.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES_DIR = Path("results")
FIG_DIR = RES_DIR / "figures"
GENES = ["TP53", "KRAS", "PIK3CA"]
K_VALUES = [5, 10]
N_PERM = 10_000
PLDDT_MIN = 70
CONTACT_A = 5.0
SEED = 42


def mean_pairwise(xyz):
    d = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(len(xyz), k=1)
    return d[iu].mean()


def load(gene):
    df = pd.read_csv(RES_DIR / f"{gene}_structure_residues.csv")
    if df["plddt"].notna().any():
        before = len(df)
        df = df[df["plddt"] >= PLDDT_MIN]
        print(f"  {gene}: kept {len(df)}/{before} residues with pLDDT >= {PLDDT_MIN}")
    return df


def clustering_test():
    rng = np.random.default_rng(SEED)
    rows, null_store = [], {}
    for gene in GENES:
        df = load(gene)
        xyz_all = df[["ca_x", "ca_y", "ca_z"]].to_numpy()
        for k in K_VALUES:
            hot = df[df["n_patients"] > 0].nlargest(k, "n_patients")
            if len(hot) < k:
                continue
            observed = mean_pairwise(hot[["ca_x", "ca_y", "ca_z"]].to_numpy())
            null = np.array([mean_pairwise(xyz_all[rng.choice(len(df), k, replace=False)])
                             for _ in range(N_PERM)])
            p = (np.sum(null <= observed) + 1) / (N_PERM + 1)
            rows.append({"gene": gene, "k_hotspots": k,
                         "hotspots": ", ".join(str(x) for x in hot["position"]),
                         "observed_mean_dist_A": round(observed, 2),
                         "random_mean_dist_A": round(null.mean(), 2),
                         "ratio_obs_vs_random": round(observed / null.mean(), 2),
                         "p_value": p})
            null_store[(gene, k)] = (null, observed, p)
    out = pd.DataFrame(rows)
    out.to_csv(RES_DIR / "clustering_test.csv", index=False)
    print("\nTest 1 - 3D clustering of hotspots")
    print(out[["gene", "k_hotspots", "observed_mean_dist_A", "random_mean_dist_A",
               "ratio_obs_vs_random", "p_value"]].to_string(index=False))

    # figure: K = 10 null distributions
    fig, axes = plt.subplots(1, len(GENES), figsize=(12, 3.6))
    for ax, gene in zip(axes, GENES):
        key = (gene, 10)
        if key not in null_store:
            ax.set_visible(False)
            continue
        null, obs, p = null_store[key]
        ax.hist(null, bins=50, color="#bdbdbd")
        ax.axvline(obs, color="#c0392b", lw=2)
        ax.text(obs, ax.get_ylim()[1] * 0.92, f"  observed\n  p = {p:.4f}",
                color="#c0392b", fontsize=8, va="top")
        ax.set_title(gene, loc="left", fontweight="bold")
        ax.set_xlabel("Mean pairwise CA distance (A)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Random residue sets")
    fig.suptitle("Are the top 10 hotspots closer together in 3D than random residues?",
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "clustering_test.png", dpi=300)
    plt.close(fig)


def distance_test(top_n=10):
    rows = []
    for gene in GENES:
        df = pd.read_csv(RES_DIR / f"{gene}_structure_residues.csv")
        if df["min_dist_to_partner"].isna().all():
            continue
        partner = df["partner"].iloc[0]
        hot = df[df["n_patients"] > 0].nlargest(top_n, "n_patients").copy()
        hot["label"] = hot["top_change"].str[0] + hot["position"].astype(str)
        hot["contact"] = hot["min_dist_to_partner"] <= CONTACT_A
        all_frac = (df["min_dist_to_partner"] <= CONTACT_A).mean()
        print(f"\n{gene}: distance of top {top_n} hotspots to {partner}")
        print(f"  {hot['contact'].sum()}/{len(hot)} hotspots within {CONTACT_A} A "
              f"(vs {100 * all_frac:.0f}% of all residues in the structure)")
        for _, r in hot.iterrows():
            print(f"    {r['label']:<6} n={r['n_patients']:<4} {r['min_dist_to_partner']:5.1f} A"
                  f"  {'CONTACT' if r['contact'] else ''}")
        hot["partner"] = partner
        hot["pct_all_residues_in_contact"] = round(100 * all_frac, 1)
        rows.append(hot)
    if not rows:
        return
    out = pd.concat(rows)
    out[["gene", "label", "position", "n_patients", "partner", "min_dist_to_partner",
         "contact", "pct_all_residues_in_contact"]].to_csv(
        RES_DIR / "distance_to_partner.csv", index=False)

    fig, axes = plt.subplots(1, len(rows), figsize=(6 * len(rows), 3.8))
    axes = np.atleast_1d(axes)
    for ax, hot in zip(axes, rows):
        hot = hot.sort_values("n_patients", ascending=False)
        colors = ["#c0392b" if c else "#7f8c8d" for c in hot["contact"]]
        ax.bar(hot["label"], hot["min_dist_to_partner"], color=colors)
        ax.axhline(CONTACT_A, ls="--", color="black", lw=0.8)
        ax.text(len(hot) - 0.5, CONTACT_A, f" {CONTACT_A:g} A contact", va="bottom",
                ha="right", fontsize=8)
        ax.set_title(f"{hot['gene'].iloc[0]}: distance to {hot['partner'].iloc[0]}",
                     loc="left", fontweight="bold")
        ax.set_ylabel("Minimum distance (A)")
        ax.tick_params(axis="x", rotation=45)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Top hotspots ordered by number of patients (red = direct contact)",
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "distance_to_partner.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    clustering_test()
    distance_test()
    print("\nDone. See results/clustering_test.csv, results/distance_to_partner.csv "
          "and results/figures/")
