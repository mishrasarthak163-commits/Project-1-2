"""
07_domain_finder.py
───────────────────
Automatically detects structured domain boundaries in AlphaFold
models using pLDDT score profiles.

The idea is simple: a stretch of residues where pLDDT stays above
a threshold (default: 70) for at least a minimum length (default: 30
residues) is likely a folded domain. Regions where pLDDT dips below
the threshold are likely flexible linkers between domains.

This is a computational alternative to running PFAM/InterPro domain
searches — it uses structural confidence rather than sequence similarity
to detect domain architecture.

Usage:
  python 07_domain_finder.py

Requirements:
  Run 05_fetch_alphafold.py and 06_extract_plddt.py first.

Output:
  Printed domain table for each protein
  ../results/domains_<protein>.csv
"""

import os
import csv
import numpy as np
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

DOMAIN_THRESHOLD  = 70   # pLDDT must be above this to count as "structured"
MIN_DOMAIN_LENGTH = 30   # minimum residues to qualify as a domain

# ── Core logic ────────────────────────────────────────────────────────────────

def extract_plddt_simple(pdb_path: str):
    """Extracts (residue_numbers, plddt_scores) from an AlphaFold PDB."""
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model     = structure[0]

    res_nums, plddt = [], []
    for chain in model:
        for residue in chain:
            if residue.get_id()[0] == " " and "CA" in residue:
                res_nums.append(residue.get_id()[1])
                plddt.append(residue["CA"].get_bfactor())

    return np.array(res_nums), np.array(plddt)


def find_domains(res_nums: np.ndarray, plddt: np.ndarray,
                 threshold: float = DOMAIN_THRESHOLD,
                 min_length: int  = MIN_DOMAIN_LENGTH) -> list:
    """
    Scans the pLDDT profile to find contiguous structured regions.

    Algorithm:
      1. Walk along the sequence residue by residue.
      2. When pLDDT crosses above `threshold`, start a new domain.
      3. When it drops below `threshold`, end the domain.
      4. Keep the domain only if it has >= `min_length` residues.

    This is similar to the approach used by EBI's SIFTS database for
    mapping structural domains to sequences.

    Args:
        res_nums:   array of residue sequence numbers
        plddt:      array of pLDDT scores (same length)
        threshold:  minimum pLDDT to qualify as structured
        min_length: minimum consecutive residues above threshold

    Returns:
        List of dicts, each describing one domain:
          {domain_n, start, end, length, mean_plddt, min_plddt, max_plddt}
    """
    domains    = []
    in_domain  = False
    start_idx  = None

    for i, (res, score) in enumerate(zip(res_nums, plddt)):

        if score >= threshold and not in_domain:
            # Entering a structured region
            in_domain = True
            start_idx = i

        elif score < threshold and in_domain:
            # Leaving the structured region
            in_domain = False
            length = i - start_idx

            if length >= min_length:
                region_plddt = plddt[start_idx:i]
                domains.append({
                    "domain_n":  len(domains) + 1,
                    "start":     int(res_nums[start_idx]),
                    "end":       int(res_nums[i - 1]),
                    "length":    length,
                    "mean_plddt":round(float(region_plddt.mean()), 1),
                    "min_plddt": round(float(region_plddt.min()), 1),
                    "max_plddt": round(float(region_plddt.max()), 1),
                })

    # Handle domain that extends to the very end of the sequence
    if in_domain:
        length = len(res_nums) - start_idx
        if length >= min_length:
            region_plddt = plddt[start_idx:]
            domains.append({
                "domain_n":  len(domains) + 1,
                "start":     int(res_nums[start_idx]),
                "end":       int(res_nums[-1]),
                "length":    length,
                "mean_plddt":round(float(region_plddt.mean()), 1),
                "min_plddt": round(float(region_plddt.min()), 1),
                "max_plddt": round(float(region_plddt.max()), 1),
            })

    return domains


def find_disordered_regions(res_nums: np.ndarray, plddt: np.ndarray,
                             threshold: float = 50,
                             min_length: int  = 10) -> list:
    """
    Finds regions with pLDDT < threshold — likely intrinsically disordered.

    These are the biologically interesting regions that:
    - Cannot be studied by X-ray crystallography (no fixed shape)
    - May be protease cleavage sites (e.g. LOXL2 propeptide)
    - May be protein-protein interaction interfaces in their free state
    - Are increasingly recognised as potential drug targets
    """
    disordered = []
    in_region  = False
    start_idx  = None

    for i, (res, score) in enumerate(zip(res_nums, plddt)):
        if score < threshold and not in_region:
            in_region = True
            start_idx = i
        elif score >= threshold and in_region:
            in_region = False
            length = i - start_idx
            if length >= min_length:
                disordered.append({
                    "start":  int(res_nums[start_idx]),
                    "end":    int(res_nums[i - 1]),
                    "length": length,
                    "mean_plddt": round(float(plddt[start_idx:i].mean()), 1),
                })

    if in_region:
        length = len(res_nums) - start_idx
        if length >= min_length:
            disordered.append({
                "start":  int(res_nums[start_idx]),
                "end":    int(res_nums[-1]),
                "length": length,
                "mean_plddt": round(float(plddt[start_idx:].mean()), 1),
            })

    return disordered


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    for name in PROTEINS:
        pdb_path = os.path.join(DATA_DIR, f"{name}_AF.pdb")

        if not os.path.exists(pdb_path):
            print(f"{name}: file missing — run 05_fetch_alphafold.py first\n")
            continue

        res_nums, plddt = extract_plddt_simple(pdb_path)
        domains         = find_domains(res_nums, plddt)
        disordered      = find_disordered_regions(res_nums, plddt)

        # ── Print results ─────────────────────────────────────────────────────
        print(f"\n{'='*62}")
        print(f" {name}  ({len(plddt)} residues total)")
        print(f"{'='*62}")

        print(f"\n  Structured domains (pLDDT > {DOMAIN_THRESHOLD}, "
              f"≥ {MIN_DOMAIN_LENGTH} residues):\n")
        print(f"  {'Domain':<8} {'Start':<7} {'End':<7} "
              f"{'Length':<8} {'Mean pLDDT':<13} {'Range'}")
        print(f"  {'-'*55}")

        for d in domains:
            plddt_range = f"{d['min_plddt']}–{d['max_plddt']}"
            print(f"  {d['domain_n']:<8} {d['start']:<7} {d['end']:<7} "
                  f"{d['length']:<8} {d['mean_plddt']:<13} {plddt_range}")

        print(f"\n  Disordered regions (pLDDT < 50, ≥ 10 residues):\n")
        if disordered:
            for r in disordered:
                print(f"  Residues {r['start']}–{r['end']}  "
                      f"({r['length']} aa, mean pLDDT: {r['mean_plddt']})")
        else:
            print("  None detected above minimum length threshold")

        # ── Save CSV ──────────────────────────────────────────────────────────
        csv_path = os.path.join(RESULTS_DIR, f"domains_{name}.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=domains[0].keys() if domains else [])
            if domains:
                writer.writeheader()
                writer.writerows(domains)
        print(f"\n  Saved → {csv_path}")

    print("\nDone. Use these domain boundaries in PyMOL to select and colour domains.")


if __name__ == "__main__":
    main()
