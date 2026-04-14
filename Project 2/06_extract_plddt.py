"""
06_extract_plddt.py
───────────────────
Extracts per-residue pLDDT confidence scores from AlphaFold PDB files
and generates a summary report and confidence plot for each protein.

What is pLDDT?
  pLDDT = predicted Local Distance Difference Test.
  It is AlphaFold's own measure of how confident it is about each
  residue's predicted position. Stored in the B-factor column of the
  PDB file (a repurposing of a field normally used for atomic vibration).

  Score interpretation:
    > 90  : Very high — trust atomic positions; as good as X-ray
    70–90 : Confident — backbone correct; minor side-chain uncertainty
    50–70 : Low — general shape likely right; flexible region
    < 50  : Very low — likely intrinsically disordered (no fixed shape)

  Importantly: a low pLDDT does NOT mean AlphaFold failed. It means
  the protein itself is flexible/disordered there. This is biologically
  meaningful — disordered regions are often cleavage sites, binding
  interfaces, or regulatory regions.

Usage:
  python 06_extract_plddt.py

Requirements:
  Run 05_fetch_alphafold.py first.

Output:
  ../results/plddt_summary.csv
  ../results/plddt_<protein>.png  (one plot per protein)
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import PDBParser

# ── Config ────────────────────────────────────────────────────────────────────

DATA_DIR    = "../data/alphafold"
RESULTS_DIR = "../results"

PROTEINS = {
    "LOXL2":  "Q9UBU1",
    "CTHRC1": "Q96CG8",
    "ITGA11": "Q9Y6M5",
    "SFRP2":  "Q96HF1",
    "PDGFRB": "P09619",
}

# Official AlphaFold colour thresholds
PLDDT_COLOURS = [
    (90,  "#006EFF", "Very high  (>90)"),
    (70,  "#65CBF3", "Confident  (70–90)"),
    (50,  "#FFDB13", "Low        (50–70)"),
    (0,   "#FF7D45", "Very low   (<50)"),
]

# ── Core extraction ───────────────────────────────────────────────────────────

def extract_plddt(pdb_path: str, protein_name: str):
    """
    Parses an AlphaFold PDB file and returns per-residue pLDDT scores.

    AlphaFold stores pLDDT in the B-factor column of each ATOM record.
    We extract the value from the CA (alpha-carbon) atom of each residue,
    giving us exactly one score per residue.

    Returns:
        (residue_numbers, plddt_scores) as numpy arrays
    """
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure(protein_name, pdb_path)
    model     = structure[0]

    res_nums = []
    plddt    = []

    for chain in model:
        for residue in chain:
            # Skip HETATM records (water, ligands)
            if residue.get_id()[0] != " ":
                continue
            # Use the CA atom as representative for the residue
            if "CA" in residue:
                res_nums.append(residue.get_id()[1])
                plddt.append(residue["CA"].get_bfactor())

    return np.array(res_nums), np.array(plddt)


def summarise_plddt(name: str, res_nums: np.ndarray, plddt: np.ndarray) -> dict:
    """Returns a dict of summary statistics for one protein."""
    n = len(plddt)
    return {
        "protein":           name,
        "total_residues":    n,
        "mean_plddt":        round(float(plddt.mean()), 1),
        "median_plddt":      round(float(np.median(plddt)), 1),
        "pct_very_high":     round(float((plddt > 90).sum() / n * 100), 1),
        "pct_confident":     round(float(((plddt > 70) & (plddt <= 90)).sum() / n * 100), 1),
        "pct_low":           round(float(((plddt > 50) & (plddt <= 70)).sum() / n * 100), 1),
        "pct_disordered":    round(float((plddt <= 50).sum() / n * 100), 1),
        "n_disordered":      int((plddt < 50).sum()),
    }


# ── Plotting ──────────────────────────────────────────────────────────────────

def colour_by_plddt(score: float) -> str:
    """Returns the AlphaFold colour for a given pLDDT score."""
    for threshold, colour, _ in PLDDT_COLOURS:
        if score >= threshold:
            return colour
    return PLDDT_COLOURS[-1][1]


def plot_plddt(name: str, res_nums: np.ndarray, plddt: np.ndarray,
               out_path: str):
    """
    Plots a per-residue pLDDT scatter plot coloured by confidence level.
    This is the same style as the official AlphaFold database plots.
    """
    colours = [colour_by_plddt(s) for s in plddt]

    fig, ax = plt.subplots(figsize=(12, 4))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")

    ax.scatter(res_nums, plddt, c=colours, s=2.5, alpha=0.9, linewidths=0)

    # Threshold lines
    for thresh, col, label in [(90, "#006EFF", ">90"), (70, "#65CBF3", "70"),
                                (50, "#FFDB13", "50")]:
        ax.axhline(y=thresh, color=col, linestyle="--", linewidth=0.8, alpha=0.6)

    ax.set_xlabel("Residue number", fontsize=10)
    ax.set_ylabel("pLDDT score", fontsize=10)
    ax.set_title(f"AlphaFold2 pLDDT — {name} ({len(plddt)} residues)",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(0, 105)
    ax.set_xlim(res_nums[0] - 5, res_nums[-1] + 5)

    # Legend
    legend_patches = [mpatches.Patch(color=col, label=lbl)
                      for _, col, lbl in PLDDT_COLOURS]
    ax.legend(handles=legend_patches, fontsize=8, loc="lower right",
              title="Confidence", title_fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    all_summaries = []

    print("Extracting pLDDT scores...\n")
    print(f"  {'Protein':<10} {'Residues':<10} {'Mean pLDDT':<13} "
          f"{'% >70':<10} {'% <50 (disordered)'}")
    print(f"  {'-'*60}")

    for name in PROTEINS:
        pdb_path = os.path.join(DATA_DIR, f"{name}_AF.pdb")

        if not os.path.exists(pdb_path):
            print(f"  {name:<10} missing — run 05_fetch_alphafold.py first")
            continue

        res_nums, plddt = extract_plddt(pdb_path, name)
        summary         = summarise_plddt(name, res_nums, plddt)
        all_summaries.append(summary)

        confident_pct   = summary["pct_very_high"] + summary["pct_confident"]
        print(f"  {name:<10} {summary['total_residues']:<10} "
              f"{summary['mean_plddt']:<13} "
              f"{confident_pct:<10.1f} "
              f"{summary['pct_disordered']}%")

        # Plot
        plot_path = os.path.join(RESULTS_DIR, f"plddt_{name}.png")
        plot_plddt(name, res_nums, plddt, plot_path)

    # Save CSV summary
    if all_summaries:
        csv_path = os.path.join(RESULTS_DIR, "plddt_summary.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_summaries[0].keys())
            writer.writeheader()
            writer.writerows(all_summaries)
        print(f"\nSummary CSV saved → {csv_path}")
        print(f"pLDDT plots saved → {RESULTS_DIR}/plddt_<protein>.png")


if __name__ == "__main__":
    main()
