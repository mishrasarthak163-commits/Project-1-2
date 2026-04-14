"""
08_superimpose_structures.py
────────────────────────────
Superimposes the AlphaFold2 model of LOXL2 onto its partial
experimental crystal structure (PDB: 5ZE3) and computes the RMSD.

What is superimposition?
  Think of it like overlaying two photographs of the same building —
  one taken by a professional (experimental structure) and one drawn
  by an AI (AlphaFold). By rotating and translating one image to best
  match the other, we can measure how similar they are.

  RMSD (Root Mean Square Deviation) is the distance metric:
    < 1.0 Å  = near-identical
    1–2 Å    = excellent agreement
    2–4 Å    = good, minor differences in loops
    > 4 Å    = significant divergence

  We compare only the residues present in BOTH structures (the
  experimentally resolved region). AlphaFold gives extra coverage
  for the propeptide and signal peptide regions not in the crystal.

Usage:
  python 08_superimpose_structures.py

Requirements:
  - Run 03_fetch_pdb_structures.py (for the experimental PDB 5ZE3)
  - Run 05_fetch_alphafold.py (for LOXL2_AF.pdb)

Output:
  LOXL2_AF_superimposed.pdb  (AlphaFold model in experimental frame)
  Printed RMSD and alignment statistics
"""

import os
from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.PDB.PDBList import PDBList

# ── Config ────────────────────────────────────────────────────────────────────

DATA_DIR_P1  = "../project1_cirrhosis_proteins/data"   # has experimental PDBs
DATA_DIR_AF  = "../data/alphafold"                     # has AlphaFold models
RESULTS_DIR  = "../results"

EXP_PDB_ID   = "5ZE3"   # partial LOXL2 crystal structure (LOX domain only)
AF_PDB_NAME  = "LOXL2"

# ── Helper functions ──────────────────────────────────────────────────────────

def get_ca_atoms(structure, chain_id: str = "A") -> dict:
    """
    Returns a dict mapping residue_number → CA_atom for one chain.
    Only includes standard amino acid residues.
    """
    model   = structure[0]
    chain   = model[chain_id]
    ca_dict = {}
    for residue in chain:
        if residue.get_id()[0] == " " and "CA" in residue:
            ca_dict[residue.get_id()[1]] = residue["CA"]
    return ca_dict


def find_common_residues(ca_exp: dict, ca_af: dict) -> list:
    """
    Returns sorted list of residue numbers present in both structures.
    These are the only residues we can fairly compare.
    """
    common = sorted(set(ca_exp.keys()) & set(ca_af.keys()))
    return common


def compute_rmsd_report(exp_structure, af_structure) -> dict:
    """
    Runs structural superimposition and returns statistics.

    BioPython's Superimposer:
    1. Finds the optimal rotation + translation to minimise RMSD
    2. Applies that transformation to the moving structure (AlphaFold)
    3. Returns the RMSD of the optimal fit
    """
    # Get CA atoms from chain A of both structures
    ca_exp = get_ca_atoms(exp_structure, "A")
    ca_af  = get_ca_atoms(af_structure,  "A")

    common = find_common_residues(ca_exp, ca_af)

    if len(common) < 10:
        return {"error": f"Only {len(common)} common residues — too few to superimpose"}

    fixed_atoms  = [ca_exp[i] for i in common]
    moving_atoms = [ca_af[i]  for i in common]

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    sup.apply(af_structure.get_atoms())   # rotate entire AF model

    return {
        "n_common_residues": len(common),
        "first_residue":     common[0],
        "last_residue":      common[-1],
        "rmsd_angstroms":    round(sup.rms, 3),
        "rotation_matrix":   sup.rotran[0].tolist(),   # 3×3 rotation
        "translation_vec":   sup.rotran[1].tolist(),   # 3-element vector
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    parser = PDBParser(QUIET=True)

    # ── Download experimental structure if needed ─────────────────────────────
    exp_pdb_path = os.path.join(DATA_DIR_P1, f"pdb{EXP_PDB_ID.lower()}.ent")
    if not os.path.exists(exp_pdb_path):
        print(f"Downloading experimental structure {EXP_PDB_ID}...")
        os.makedirs(DATA_DIR_P1, exist_ok=True)
        pdbl = PDBList(verbose=False)
        pdbl.retrieve_pdb_file(EXP_PDB_ID, pdir=DATA_DIR_P1, file_type="pdb")

    # ── Load structures ───────────────────────────────────────────────────────
    af_pdb_path = os.path.join(DATA_DIR_AF, f"{AF_PDB_NAME}_AF.pdb")

    for path, label in [(exp_pdb_path, "experimental"), (af_pdb_path, "AlphaFold")]:
        if not os.path.exists(path):
            print(f"Missing {label} structure: {path}")
            print("Please run the prerequisite scripts first.")
            return

    print(f"\nLoading structures...")
    exp_structure = parser.get_structure("LOXL2_exp", exp_pdb_path)
    af_structure  = parser.get_structure("LOXL2_AF",  af_pdb_path)

    print(f"  Experimental (5ZE3) loaded")
    print(f"  AlphaFold (Q9UBU1) loaded")

    # ── Superimpose ───────────────────────────────────────────────────────────
    print("\nRunning superimposition (CA atoms only)...")
    result = compute_rmsd_report(exp_structure, af_structure)

    if "error" in result:
        print(f"Superimposition failed: {result['error']}")
        return

    # ── Print report ──────────────────────────────────────────────────────────
    print("\n" + "=" * 55)
    print("SUPERIMPOSITION RESULTS: LOXL2 AlphaFold vs 5ZE3")
    print("=" * 55)
    print(f"  Common residues aligned : {result['n_common_residues']}")
    print(f"  Residue range           : {result['first_residue']}–{result['last_residue']}")
    print(f"  RMSD                    : {result['rmsd_angstroms']} Å")
    print()

    rmsd = result["rmsd_angstroms"]
    if rmsd < 1.0:
        quality = "Near-identical — excellent prediction"
    elif rmsd < 2.0:
        quality = "Excellent agreement"
    elif rmsd < 4.0:
        quality = "Good — minor loop differences"
    else:
        quality = "Significant divergence in some regions"

    print(f"  Quality assessment      : {quality}")
    print()
    print("  Interpretation:")
    print("  The AlphaFold model was rotated and translated to best")
    print("  overlap with the experimental structure. The RMSD tells")
    print("  us the average distance between matching CA atoms after")
    print("  this optimal fit.")
    print("=" * 55)

    # ── Save superimposed structure ───────────────────────────────────────────
    out_path = os.path.join(RESULTS_DIR, "LOXL2_AF_superimposed.pdb")
    io = PDBIO()
    io.set_structure(af_structure)
    io.save(out_path)
    print(f"\nSuperimposed AlphaFold model saved → {out_path}")
    print("Load both pdb5ze3.ent and LOXL2_AF_superimposed.pdb in PyMOL")
    print("and run: colour marine, LOXL2_exp; colour orange, LOXL2_AF")


if __name__ == "__main__":
    main()
