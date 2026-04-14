"""
02_protparam_analysis.py
────────────────────────
Calculates physico-chemical properties for all four cirrhosis-linked
proteins using BioPython's ProteinAnalysis module.

Think of this script as running a "blood test" on each protein —
it tells us the protein's weight, charge, hydrophobicity, and how
much of it folds into helices vs. sheets.

Properties calculated:
  • Molecular weight (Da)
  • Isoelectric point (pI) — the pH at which the protein has zero net charge
  • Instability index — < 40 means the protein is stable in solution
  • GRAVY score — negative = hydrophilic (water-loving), positive = hydrophobic
  • Secondary structure fractions — alpha helix, beta sheet, random coil

Usage:
  python 02_protparam_analysis.py

Requirements:
  Run 01_fetch_sequences.py first to generate the FASTA files.

Output:
  Printed summary table + ../results/protparam_summary.csv
"""

import os
import csv
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ── Config ────────────────────────────────────────────────────────────────────

DATA_DIR    = "../data"
RESULTS_DIR = "../results"

PROTEINS = ["TGF_beta1", "ACTA2", "COL1A1", "TIMP1"]

# ── Analysis ──────────────────────────────────────────────────────────────────

def analyse_protein(fasta_path: str) -> dict:
    """
    Reads a FASTA file and computes key physico-chemical properties.

    BioPython's ProteinAnalysis class does all the heavy lifting here.
    It implements the same algorithms as the ExPASy ProtParam web server,
    but entirely in Python so we can automate it for many proteins at once.

    Args:
        fasta_path: path to a .fasta file containing one protein sequence

    Returns:
        dict of computed properties
    """
    record   = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq)

    # Remove any ambiguous amino acid characters (X, B, Z) — ProteinAnalysis
    # will raise an error if they are present
    clean_seq = "".join(aa for aa in sequence if aa not in "XBZ*")

    pa = ProteinAnalysis(clean_seq)

    mw          = pa.molecular_weight()
    pi          = pa.isoelectric_point()
    instability = pa.instability_index()
    gravy       = pa.gravy()

    # secondary_structure_fraction() returns (helix, turn, sheet) fractions
    helix, turn, sheet = pa.secondary_structure_fraction()

    # Amino acid composition as percentages
    aa_comp = pa.get_amino_acids_percent()

    # Top 3 most abundant amino acids
    top3 = sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)[:3]
    top3_str = ", ".join(f"{aa}:{pct*100:.1f}%" for aa, pct in top3)

    return {
        "sequence_length":  len(clean_seq),
        "molecular_weight": round(mw, 1),
        "pi":               round(pi, 2),
        "instability_index":round(instability, 2),
        "stability":        "Stable" if instability < 40 else "Unstable",
        "gravy":            round(gravy, 3),
        "hydro_character":  "Hydrophilic" if gravy < 0 else "Hydrophobic",
        "helix_pct":        round(helix * 100, 1),
        "sheet_pct":        round(sheet * 100, 1),
        "coil_pct":         round((1 - helix - turn - sheet) * 100, 1),
        "top3_aa":          top3_str,
    }


def print_summary(results: dict):
    """Prints a formatted summary table to the terminal."""

    props = [
        ("Sequence length (aa)",  "sequence_length"),
        ("Molecular weight (Da)", "molecular_weight"),
        ("Isoelectric point (pI)","pi"),
        ("Instability index",     "instability_index"),
        ("Stability",             "stability"),
        ("GRAVY score",           "gravy"),
        ("Hydro. character",      "hydro_character"),
        ("Alpha helix %",         "helix_pct"),
        ("Beta sheet %",          "sheet_pct"),
        ("Random coil %",         "coil_pct"),
        ("Top 3 amino acids",     "top3_aa"),
    ]

    names = list(results.keys())
    col_w = 24

    # Header
    print("\n" + "=" * (col_w + len(names) * 20))
    print(f"{'Property':<{col_w}}" + "".join(f"{n:>20}" for n in names))
    print("-" * (col_w + len(names) * 20))

    for label, key in props:
        row = f"{label:<{col_w}}"
        for name in names:
            val = str(results[name][key])
            row += f"{val:>20}"
        print(row)

    print("=" * (col_w + len(names) * 20))


def save_csv(results: dict, out_path: str):
    """Saves results as a CSV file for further analysis."""
    all_keys = list(next(iter(results.values())).keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["property"] + list(results.keys()))
        for key in all_keys:
            row = [key] + [results[name][key] for name in results]
            writer.writerow(row)
    print(f"\nCSV saved → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    results = {}

    print("Running ProteinAnalysis on all four proteins...\n")

    for name in PROTEINS:
        fasta_path = os.path.join(DATA_DIR, f"{name}.fasta")

        if not os.path.exists(fasta_path):
            print(f"  MISSING: {fasta_path} — run 01_fetch_sequences.py first")
            continue

        print(f"  Analysing {name}...", end=" ", flush=True)
        results[name] = analyse_protein(fasta_path)
        print("done")

    print_summary(results)

    csv_path = os.path.join(RESULTS_DIR, "protparam_summary.csv")
    save_csv(results, csv_path)


if __name__ == "__main__":
    main()
