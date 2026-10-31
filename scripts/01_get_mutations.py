"""
01_get_mutations.py
-------------------
Download missense mutations for TP53, KRAS and PIK3CA from all
TCGA PanCancer Atlas studies on cBioPortal, then count how many
patients carry a mutation at each protein position.

Outputs:
  data/raw/<GENE>_mutations_raw.csv          all mutations (every type)
  data/processed/<GENE>_missense_counts.csv  patients per residue (missense only)

Run from the project folder:
  python scripts/01_get_mutations.py
"""

import time
from pathlib import Path

import pandas as pd
import requests

API = "https://www.cbioportal.org/api"

# Gene symbol -> Entrez Gene ID (what cBioPortal uses)
GENES = {"TP53": 7157, "KRAS": 3845, "PIK3CA": 5290}

RAW_DIR = Path("data/raw")
PROC_DIR = Path("data/processed")
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR.mkdir(parents=True, exist_ok=True)


def get_json(url, params=None, tries=3):
    """GET a URL and return JSON, retrying a few times if the server hiccups."""
    for attempt in range(1, tries + 1):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            if attempt == tries:
                print(f"   ! failed after {tries} tries: {e}")
                return None
            time.sleep(2 * attempt)


def get_pancan_studies():
    """Return the IDs of all TCGA PanCancer Atlas studies (one per cancer type)."""
    studies = get_json(f"{API}/studies", params={"projection": "SUMMARY"})
    ids = sorted(s["studyId"] for s in studies
                 if s["studyId"].endswith("_tcga_pan_can_atlas_2018"))
    print(f"Found {len(ids)} TCGA PanCancer Atlas studies")
    return ids


def get_mutations(study_id, entrez_id):
    """All mutations in one gene for one study."""
    url = f"{API}/molecular-profiles/{study_id}_mutations/mutations"
    params = {
        "sampleListId": f"{study_id}_all",
        "entrezGeneId": entrez_id,
        "projection": "DETAILED",
    }
    data = get_json(url, params)
    return data or []


def main():
    studies = get_pancan_studies()

    for gene, entrez in GENES.items():
        print(f"\n=== {gene} ===")
        rows = []
        for study in studies:
            muts = get_mutations(study, entrez)
            for m in muts:
                rows.append({
                    "study": study,
                    "cancer_type": study.replace("_tcga_pan_can_atlas_2018", ""),
                    "patient": m.get("patientId"),
                    "sample": m.get("sampleId"),
                    "protein_change": m.get("proteinChange"),
                    "mutation_type": m.get("mutationType"),
                    "position": m.get("proteinPosStart"),
                })
            print(f"  {study:<40} {len(muts):>4} mutations")

        df = pd.DataFrame(rows)
        df.to_csv(RAW_DIR / f"{gene}_mutations_raw.csv", index=False)

        # Keep missense only, with a valid position
        mis = df[(df["mutation_type"] == "Missense_Mutation") & (df["position"] > 0)].copy()
        mis["position"] = mis["position"].astype(int)

        # Count each patient once per position (avoids double-counting multiple samples)
        mis = mis.drop_duplicates(subset=["patient", "position"])

        counts = (mis.groupby("position")
                     .agg(n_patients=("patient", "nunique"),
                          top_change=("protein_change", lambda x: x.value_counts().index[0]))
                     .reset_index()
                     .sort_values("position"))
        counts.insert(0, "gene", gene)
        counts.to_csv(PROC_DIR / f"{gene}_missense_counts.csv", index=False)

        print(f"\n  Total mutations: {len(df)} | missense patients-positions: {len(mis)}")
        print("  Top 5 hotspots:")
        top = counts.sort_values("n_patients", ascending=False).head(5)
        for _, r in top.iterrows():
            print(f"    position {r['position']:>5}  {r['top_change']:<8} {r['n_patients']} patients")

    print("\nDone. Check data/raw/ and data/processed/")


if __name__ == "__main__":
    main()
