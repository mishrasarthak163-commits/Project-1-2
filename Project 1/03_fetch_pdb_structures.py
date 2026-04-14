"""
03_fetch_pdb_structures.py
──────────────────────────
Downloads 3D crystal structures from the RCSB Protein Data Bank (PDB)
for our four cirrhosis proteins and prints a structural summary.

The PDB file format stores the x, y, z coordinates of every atom in a
protein — thousands of numbers representing a 3D sculpture at atomic
resolution. BioPython reads this and gives us a clean Python object we
can navigate like a directory tree:

  Structure
  └── Model (usually just one)
      └── Chain (A, B, C ... one per polypeptide chain)
          └── Residue (each amino acid)
              └── Atom (N, CA, CB, C, O, side-chain atoms...)

Usage:
  python 03_fetch_pdb_structures.py

Output:
  PDB files downloaded to ../data/
  Structural summary printed to terminal
"""

import os
from Bio.PDB import PDBParser
from Bio.PDB.PDBList import PDBList

# ── Config ────────────────────────────────────────────────────────────────────

DATA_DIR = "../data"

PDB_IDS = {
    "TGF_beta1": "1KTZ",   # TGF-beta1 — crystal structure of the dimer
    "ACTA2":     "3LUE",   # Actin — filament-form structure
    "COL1A1":    "1BKV",   # Collagen triple helix peptide model
    "TIMP1":     "1UEA",   # TIMP1 bound to MMP-3 (shows inhibitory mechanism)
}

# Residues that are NOT amino acids (water, ions, ligands) — skip these
NON_AA = {"HOH", "SO4", "PO4", "GOL", "EDO", "PEG", "NAG", "MAN", "FUC",
           "ZN", "CA", "MG", "CL", "NA", "K"}

# ── Functions ─────────────────────────────────────────────────────────────────

def download_structures(pdb_ids: dict, output_dir: str):
    """
    Downloads PDB files using BioPython's PDBList.

    PDBList.retrieve_pdb_file() fetches the file from RCSB and saves it
    locally. It uses the format 'pdb{id}.ent' by default.
    """
    os.makedirs(output_dir, exist_ok=True)
    pdbl = PDBList(verbose=False)

    print("Downloading PDB structures from RCSB...\n")
    for name, pdb_id in pdb_ids.items():
        print(f"  [{pdb_id}] {name}...", end=" ", flush=True)
        pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_type="pdb")
        print("done")


def summarise_structure(pdb_id: str, data_dir: str) -> dict:
    """
    Parses a PDB file and extracts key structural statistics.

    Returns a dict with:
      - number of chains
      - residue count per chain
      - total atom count
      - resolution (if available in REMARK 2)
    """
    pdb_file = os.path.join(data_dir, f"pdb{pdb_id.lower()}.ent")

    if not os.path.exists(pdb_file):
        return {"error": f"File not found: {pdb_file}"}

    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    model     = structure[0]

    chain_info    = {}
    total_atoms   = 0
    total_residues = 0

    for chain in model:
        aa_residues = [
            r for r in chain.get_residues()
            if r.get_id()[0] == " " and r.get_resname().strip() not in NON_AA
        ]
        n_atoms   = sum(len(list(r.get_atoms())) for r in aa_residues)
        chain_info[chain.id] = len(aa_residues)
        total_residues += len(aa_residues)
        total_atoms    += n_atoms

    # Also count all atoms including HETATM (ligands, water)
    all_atoms = sum(1 for _ in structure.get_atoms())

    return {
        "chains":            list(chain_info.keys()),
        "residues_per_chain": chain_info,
        "total_aa_residues": total_residues,
        "total_atoms":       all_atoms,
    }


def print_structural_report(pdb_ids: dict, data_dir: str):
    """Prints a formatted structural summary for all proteins."""

    print("\n" + "=" * 65)
    print("STRUCTURAL SUMMARY — PDB STRUCTURES")
    print("=" * 65)

    for name, pdb_id in pdb_ids.items():
        info = summarise_structure(pdb_id, data_dir)

        if "error" in info:
            print(f"\n{name} ({pdb_id}): {info['error']}")
            continue

        print(f"\n{name}  (PDB: {pdb_id})")
        print(f"  Chains          : {', '.join(info['chains'])}")

        for chain_id, n_res in info["residues_per_chain"].items():
            print(f"  Chain {chain_id}         : {n_res} amino acid residues")

        print(f"  Total residues  : {info['total_aa_residues']}")
        print(f"  Total atoms     : {info['total_atoms']}")

    print("\n" + "=" * 65)
    print("Tip: Open these .ent files in PyMOL to visualise the structures.")
    print(f"Files location: {os.path.abspath(data_dir)}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    download_structures(PDB_IDS, DATA_DIR)
    print_structural_report(PDB_IDS, DATA_DIR)


if __name__ == "__main__":
    main()
