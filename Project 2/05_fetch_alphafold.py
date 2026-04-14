"""
05_fetch_alphafold.py
─────────────────────
Downloads AlphaFold2 predicted structures for five understudied
fibrosis proteins from the EBI AlphaFold database.

Background:
  AlphaFold2 (DeepMind, 2021) can predict the 3D structure of any
  protein from its amino acid sequence alone. EBI hosts these
  predictions for the entire human proteome, freely available.

  The files are in standard PDB format, with one important twist:
  the B-factor column (normally used for atomic vibration in
  X-ray structures) stores the pLDDT confidence score (0–100).
  A pLDDT > 90 means the prediction is essentially as reliable
  as an experimental structure.

Proteins targeted:
  LOXL2   (Q9UBU1) – cross-links collagen, major fibrosis driver
  CTHRC1  (Q96CG8) – collagen triple helix repeat protein
  ITGA11  (Q9Y6M5) – integrin receptor for collagen I
  SFRP2   (Q96HF1) – Wnt pathway inhibitor
  PDGFRB  (P09619) – PDGF receptor, drives HSC proliferation

Usage:
  python 05_fetch_alphafold.py

Output:
  ../data/alphafold/LOXL2_AF.pdb
  ../data/alphafold/CTHRC1_AF.pdb
  ... (one file per protein)
"""

import urllib.request
import os
import time

# ── Config ────────────────────────────────────────────────────────────────────

DATA_DIR    = "../data/alphafold"
AF_URL_TMPL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"

TARGETS = {
    "LOXL2":  "Q9UBU1",   # 774 aa — lysyl oxidase homolog
    "CTHRC1": "Q96CG8",   # 243 aa — collagen triple helix repeat
    "ITGA11": "Q9Y6M5",   # 1188 aa — integrin alpha-11
    "SFRP2":  "Q96HF1",   # 295 aa — secreted frizzled-related protein 2
    "PDGFRB": "P09619",   # 1106 aa — PDGF receptor beta
}

# ── Main ──────────────────────────────────────────────────────────────────────

def download_alphafold_model(accession: str, out_path: str) -> bool:
    """
    Downloads a single AlphaFold model PDB file from EBI.

    Args:
        accession: UniProt accession (e.g. 'Q9UBU1')
        out_path:  where to save the .pdb file

    Returns:
        True if successful, False otherwise
    """
    url = AF_URL_TMPL.format(acc=accession)
    try:
        urllib.request.urlretrieve(url, out_path)
        return True
    except urllib.error.HTTPError as e:
        print(f"    HTTP {e.code}: {e.reason}")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def count_residues_in_pdb(pdb_path: str) -> int:
    """
    Quickly counts residues by scanning for unique ATOM CA records.
    Much faster than parsing with BioPython for a simple count.
    """
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain   = line[21]
                res_num = line[22:26].strip()
                seen.add((chain, res_num))
    return len(seen)


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    print("Downloading AlphaFold2 models from EBI...\n")
    print(f"  {'Protein':<10} {'UniProt':<10} {'Status':<12} {'Residues'}")
    print(f"  {'-'*50}")

    for name, acc in TARGETS.items():
        out_path = os.path.join(DATA_DIR, f"{name}_AF.pdb")

        if os.path.exists(out_path):
            n_res = count_residues_in_pdb(out_path)
            print(f"  {name:<10} {acc:<10} {'cached':<12} {n_res} aa")
            continue

        success = download_alphafold_model(acc, out_path)

        if success:
            n_res = count_residues_in_pdb(out_path)
            print(f"  {name:<10} {acc:<10} {'downloaded':<12} {n_res} aa")
        else:
            print(f"  {name:<10} {acc:<10} {'FAILED':<12} —")

        time.sleep(0.5)   # polite pause between requests

    print(f"\nAll files saved to: {os.path.abspath(DATA_DIR)}")
    print("Next step: run 06_extract_plddt.py")


if __name__ == "__main__":
    main()
