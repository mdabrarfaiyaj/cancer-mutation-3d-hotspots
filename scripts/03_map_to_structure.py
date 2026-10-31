"""
03_map_to_structure.py
----------------------
Put the mutation counts from script 01 onto 3D protein structures.

Structures used:
  TP53   PDB 1TSR  - p53 DNA-binding domain bound to DNA (the chain that makes
                     the most DNA contacts is picked automatically)
  KRAS   PDB 4OBE  - wild-type KRAS bound to GDP + Mg (chain A)
  PIK3CA AlphaFold model AF-P42336 (no experimental full-length structure needed)

For each gene this script:
  1. downloads the structure if it is missing
  2. writes a per-residue table: mutation count, CA coordinates, AlphaFold
     confidence (pLDDT), and minimum distance to DNA / GDP
  3. writes a PDB file with the mutation count stored in the B-factor column
     (open it in ChimeraX/PyMOL and colour by B-factor)
  4. renders a PNG image with PyMOL (skipped if PyMOL is not installed)
  5. writes an interactive 3D HTML page (3Dmol.js) you can put on your website

Outputs:
  results/<GENE>_structure_residues.csv
  results/structures/<GENE>_mutcount.pdb
  results/figures/<GENE>_3D.png
  results/interactive/<GENE>_3D.html

Run from the project folder:
  python scripts/03_map_to_structure.py
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBIO, NeighborSearch, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from matplotlib import colormaps
from matplotlib.colors import to_hex

PROC_DIR = Path("data/processed")
STRUCT_DIR = Path("structures")
OUT_STRUCT = Path("results/structures")
OUT_FIG = Path("results/figures")
OUT_HTML = Path("results/interactive")
for d in (STRUCT_DIR, OUT_STRUCT, OUT_FIG, OUT_HTML):
    d.mkdir(parents=True, exist_ok=True)

DNA_NAMES = {"DA", "DC", "DG", "DT"}

GENES = {
    "TP53": {"file": "1TSR.pdb", "source": "pdb", "id": "1TSR",
             "chain": "auto_dna", "partner": "DNA"},
    "KRAS": {"file": "4OBE.pdb", "source": "pdb", "id": "4OBE",
             "chain": "A", "partner": "GDP"},
    "PIK3CA": {"file": "AF-P42336.pdb", "source": "alphafold", "id": "P42336",
               "chain": "A", "partner": None},
}

N_LABELS = 6  # hotspots labelled in images / html


# ---------------------------------------------------------------- download
def download_structure(cfg):
    path = STRUCT_DIR / cfg["file"]
    if path.exists():
        return path
    if cfg["source"] == "pdb":
        url = f"https://files.rcsb.org/download/{cfg['id']}.pdb"
    else:
        meta = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{cfg['id']}",
                            timeout=60).json()
        url = meta[0]["pdbUrl"]
    print(f"  downloading {url}")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    path.write_text(r.text)
    return path


# ---------------------------------------------------------------- helpers
def protein_residues(chain):
    return [r for r in chain if r.id[0] == " " and is_aa(r, standard=True)]


def partner_atoms(model, chain_id, partner):
    """Atoms of the functional partner: DNA (any chain) or GDP (in the chosen chain)."""
    atoms = []
    for chain in model:
        for res in chain:
            name = res.get_resname().strip()
            if partner == "DNA" and name in DNA_NAMES:
                atoms.extend(res.get_atoms())
            if partner == "GDP" and chain.id == chain_id and name == "GDP":
                atoms.extend(res.get_atoms())
    return atoms


def pick_chain(model, cfg):
    """Use the configured chain, or (for TP53) the chain with the most DNA contacts."""
    if cfg["chain"] != "auto_dna":
        return cfg["chain"]
    dna = partner_atoms(model, None, "DNA")
    ns = NeighborSearch(dna)
    best, best_n = None, -1
    for chain in model:
        res = protein_residues(chain)
        if not res:
            continue
        n = sum(1 for r in res for a in r if ns.search(a.coord, 4.0))
        print(f"    chain {chain.id}: {n} protein atoms within 4 A of DNA")
        if n > best_n:
            best, best_n = chain.id, n
    return best


class KeepSelect(Select):
    """Keep the chosen protein chain plus its DNA/GDP partner; drop water and other copies."""

    def __init__(self, chain_id, partner):
        self.chain_id, self.partner = chain_id, partner

    def accept_residue(self, res):
        name = res.get_resname().strip()
        chain_id = res.get_parent().id
        if res.id[0] == " " and is_aa(res, standard=True):
            return chain_id == self.chain_id
        if self.partner == "DNA" and name in DNA_NAMES:
            return True
        if self.partner == "GDP" and chain_id == self.chain_id and name in {"GDP", "MG"}:
            return True
        return False


# ---------------------------------------------------------------- per gene
def process_gene(gene, cfg):
    print(f"\n=== {gene} ===")
    path = download_structure(cfg)
    model = PDBParser(QUIET=True).get_structure(gene, path)[0]
    chain_id = pick_chain(model, cfg)
    chain = model[chain_id]
    print(f"  using chain {chain_id} of {path.name}")

    counts = pd.read_csv(PROC_DIR / f"{gene}_missense_counts.csv")
    count_map = dict(zip(counts["position"], counts["n_patients"]))
    change_map = dict(zip(counts["position"], counts["top_change"]))

    partner = partner_atoms(model, chain_id, cfg["partner"]) if cfg["partner"] else []
    ns = NeighborSearch(partner) if partner else None
    partner_xyz = np.array([a.coord for a in partner]) if partner else None
    if cfg["partner"] and not partner:
        print(f"  ! no {cfg['partner']} atoms found - distances will be empty")

    rows = []
    for res in protein_residues(chain):
        pos = res.id[1]
        ca = res["CA"].coord if "CA" in res else None
        if ca is None:
            continue
        n = int(count_map.get(pos, 0))
        dist = np.nan
        if partner_xyz is not None:
            res_xyz = np.array([a.coord for a in res])
            diff = res_xyz[:, None, :] - partner_xyz[None, :, :]
            dist = float(np.sqrt((diff ** 2).sum(-1)).min())
        plddt = float(res["CA"].get_bfactor()) if cfg["source"] == "alphafold" else np.nan
        rows.append({"gene": gene, "position": pos, "residue": res.get_resname(),
                     "n_patients": n, "top_change": change_map.get(pos, ""),
                     "ca_x": ca[0], "ca_y": ca[1], "ca_z": ca[2],
                     "plddt": plddt, "min_dist_to_partner": dist,
                     "partner": cfg["partner"] or ""})
        # store the mutation count in the B-factor column
        for atom in res:
            atom.set_bfactor(n)

    # zero B-factors for DNA / ligand so colouring by B shows only mutations
    for ch in model:
        for res in ch:
            if not (res.id[0] == " " and is_aa(res, standard=True)):
                for atom in res:
                    atom.set_bfactor(0)

    table = pd.DataFrame(rows)
    table.to_csv(Path("results") / f"{gene}_structure_residues.csv", index=False)

    covered = table.loc[table["n_patients"] > 0, "n_patients"].sum()
    print(f"  residues in structure: {table['position'].min()}-{table['position'].max()} "
          f"({len(table)} residues)")
    print(f"  mutations covered by this structure: {covered} of {counts['n_patients'].sum()} "
          f"({100 * covered / counts['n_patients'].sum():.0f}%)")

    out_pdb = OUT_STRUCT / f"{gene}_mutcount.pdb"
    io = PDBIO()
    io.set_structure(model)
    io.save(str(out_pdb), KeepSelect(chain_id, cfg["partner"]))

    top = table.nlargest(N_LABELS, "n_patients")
    top = top[top["n_patients"] > 0]
    render_pymol(gene, out_pdb, top, table["n_patients"].max())
    write_html(gene, out_pdb, table, top)


# ---------------------------------------------------------------- images
def residue_label(row):
    return f"{str(row['top_change'])[:1]}{row['position']}"


def render_pymol(gene, pdb_path, top, nmax):
    try:
        import pymol2
    except ImportError:
        print("  (PyMOL not installed - skipping PNG)")
        return
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(str(pdb_path), "s")
        cmd.hide("everything")
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 1)
        cmd.show("cartoon", "polymer.protein")
        # colour on a log scale so one huge hotspot (KRAS G12) doesn't wash out the rest
        cmd.alter("polymer.protein", "b=math.log10(b+1)", space={"math": math})
        cmd.spectrum("b", "white_red", "polymer.protein",
                      minimum=0, maximum=math.log10(nmax + 1))
        cmd.show("cartoon", "polymer.nucleic")
        cmd.set("cartoon_ring_mode", 3)
        cmd.color("lightorange", "polymer.nucleic")
        cmd.show("sticks", "organic")
        cmd.color("marine", "organic")
        cmd.show("spheres", "inorganic")
        for _, r in top.iterrows():
            sel = f"polymer.protein and resi {r['position']} and name CA"
            cmd.show("spheres", sel)
            cmd.set("sphere_scale", 0.9, sel)
            cmd.label(sel, f'"{residue_label(r)}"')
        cmd.set("label_size", 22)
        cmd.set("label_color", "black")
        cmd.set("label_position", (0, 2.5, 4))
        cmd.orient("all")
        cmd.zoom("all", 3)
        out = OUT_FIG / f"{gene}_3D.png"
        cmd.png(str(out), width=1800, height=1400, dpi=300, ray=1)
    print(f"  saved {out}")


# ---------------------------------------------------------------- interactive
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{gene} mutation hotspots in 3D</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 0; padding: 16px; color: #222; }}
 #viewer {{ width: 100%; height: 560px; position: relative; border: 1px solid #ddd; border-radius: 8px; }}
 .legend span {{ display: inline-block; width: 12px; height: 12px; border-radius: 50%; margin: 0 4px 0 12px; vertical-align: middle; }}
 p.note {{ font-size: 0.9em; color: #555; }}
</style></head><body>
<h2>{gene}: where cancer mutations sit in 3D</h2>
<p class="note">{note} Spheres mark mutated residues (bigger and redder = more patients,
TCGA PanCancer Atlas, 32 cancer types). Hover over a sphere for details. Drag to rotate, scroll to zoom.</p>
<div class="legend">{legend}</div>
<div id="viewer"></div>
<script>
const pdb = {pdb_json};
const muts = {muts_json};
const v = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
v.addModel(pdb, "pdb");
v.setStyle({{}}, {{cartoon: {{color: "#d9d9d9"}}}});
v.setStyle({{resn: ["DA","DC","DG","DT"]}}, {{cartoon: {{color: "#fdae6b"}}, stick: {{color: "#fdae6b", radius: 0.15}}}});
v.setStyle({{resn: ["GDP"]}}, {{stick: {{colorscheme: "blueCarbon"}}}});
v.setStyle({{resn: ["MG"]}}, {{sphere: {{color: "green", radius: 0.6}}}});
muts.forEach(m => {{
  v.addStyle({{resi: m.pos, atom: "CA"}}, {{sphere: {{color: m.color, radius: m.radius}}}});
  if (m.label) v.addLabel(m.name, {{position: {{x: m.x, y: m.y, z: m.z}}, fontSize: 13,
      fontColor: "black", backgroundColor: "white", backgroundOpacity: 0.75, inFront: true}});
}});
const byPos = Object.fromEntries(muts.map(m => [m.pos, m]));
let hover = null;
v.setHoverable({{atom: "CA"}}, true,
  function(atom) {{
    const m = byPos[atom.resi]; if (!m || hover) return;
    hover = v.addLabel(m.name + ": " + m.n + " patients (most common " + m.change + ")",
      {{position: atom, fontSize: 12, backgroundColor: "#333", fontColor: "white"}});
  }},
  function() {{ if (hover) {{ v.removeLabel(hover); hover = null; }} }});
v.zoomTo(); v.render();
</script></body></html>
"""

NOTES = {
    "TP53": "Grey = p53 DNA-binding domain (PDB 1TSR), orange = DNA.",
    "KRAS": "Grey = KRAS (PDB 4OBE), blue sticks = GDP, green = magnesium.",
    "PIK3CA": "Grey = PIK3CA AlphaFold model (AF-P42336); no partner protein shown.",
}


def write_html(gene, pdb_path, table, top):
    mut = table[table["n_patients"] > 0].copy()
    nmax = mut["n_patients"].max()
    cmap = colormaps["Reds"]
    top_pos = set(top["position"])
    muts = []
    for _, r in mut.iterrows():
        scale = math.log10(r["n_patients"] + 1) / math.log10(nmax + 1)
        muts.append({
            "pos": int(r["position"]), "n": int(r["n_patients"]),
            "name": residue_label(r), "change": str(r["top_change"]),
            "color": to_hex(cmap(0.3 + 0.7 * scale)),
            "radius": round(0.6 + 1.6 * scale, 2),
            "x": float(r["ca_x"]), "y": float(r["ca_y"]), "z": float(r["ca_z"]),
            "label": int(r["position"]) in top_pos,
        })
    legend = "".join(
        f'<span style="background:{to_hex(cmap(0.3 + 0.7 * s))}"></span>{lab}'
        for s, lab in [(0.15, "few patients"), (0.6, "more"), (1.0, "most")]
    )
    html = HTML_TEMPLATE.format(gene=gene, note=NOTES[gene], legend=legend,
                                pdb_json=json.dumps(pdb_path.read_text()),
                                muts_json=json.dumps(muts))
    out = OUT_HTML / f"{gene}_3D.html"
    out.write_text(html)
    print(f"  saved {out}")


if __name__ == "__main__":
    for gene, cfg in GENES.items():
        process_gene(gene, cfg)
    print("\nDone. Open results/interactive/*.html in your browser, "
          "and results/structures/*_mutcount.pdb in ChimeraX (colour by B-factor).")
