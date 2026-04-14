"""
04_residue_distances.py
───────────────────────
Finds residues near a key active site residue and measures
inter-residue distances within protein structures.

Why does distance matter? In a protein, two residues can only
interact (form hydrogen bonds, hydrophobic contacts, etc.) if they
are physically close — typically within 3–5 Ångströms of each other.
By measuring distances, we can:
  - Define the active site (all residues within 5 Å of the catalytic residue)
  - Identify stabilising contacts
  - Spot potential drug-binding pockets

This script analyses TIMP1 (PDB: 1UEA) — specifically how Cys1
(the inhibitory "tip") sits inside the MMP active site.

Usage:
  python 04_residue_distances.py

Requirements:
  Run 03_fetch_pdb_structures.py first.

Output:
  Printed neighbour lists + ../results/timp1_active_site.csv
"""

import os
import csv
from Bio.PDB import PDBParser, NeighborSearch

# ── Config ────────────────────────────────────────────────────────────────────

DATA_DIR    = "../data"
RESULTS_DIR = "../results"

TIMP1_PDB   = os.path.join(DATA_DIR, "pdb1uea.ent")

# ── Functions ─────────────────────────────────────────────────────────────────

def find_nearby_residues(structure, chain_id: str, res_id: int,
                         radius: float = 5.0) -> list:
    """
    Finds all residues within `radius` Ångströms of a target residue.

    Uses BioPython's NeighborSearch — an efficient spatial index that
    avoids checking every pair of atoms (which would be O(N²) slow).
    Instead it builds a KD-tree so searches run in O(log N) time.

    Args:
        structure:  a BioPython Structure object
        chain_id:   which chain to look in (e.g. 'A')
        res_id:     residue sequence number (e.g. 1 for the first residue)
        radius:     search radius in Ångströms

    Returns:
        List of nearby Residue objects (sorted by residue number)
    """
    model     = structure[0]
    all_atoms = list(model.get_atoms())
    ns        = NeighborSearch(all_atoms)

    # Get the alpha-carbon (CA) of the target residue as our search centre.
    # CA is the central atom of each amino acid backbone.
    target_res = model[chain_id][res_id]
    if "CA" not in target_res:
        raise ValueError(f"No CA atom in residue {res_id} of chain {chain_id}")

    centre = target_res["CA"].get_vector().get_array()

    nearby = ns.search(centre, radius, level="R")   # level='R' = residue level
    nearby = sorted(nearby, key=lambda r: r.get_id()[1])

    return nearby


def ca_distance(res1, res2) -> float:
    """
    Calculates the straight-line distance (Å) between the alpha-carbons
    of two residues.

    Alpha-carbons are the 'backbone' atoms — using them gives a consistent
    measure of how far apart two residues are in space, ignoring the
    variable side chains.

    Args:
        res1, res2: BioPython Residue objects

    Returns:
        Distance in Ångströms
    """
    diff = res1["CA"].get_vector() - res2["CA"].get_vector()
    return diff.norm()


def all_pairwise_distances(chain, max_residues: int = 50) -> list:
    """
    Computes CA-CA distances for all pairs of residues in a chain
    (up to max_residues to keep output manageable).

    Returns a list of (res1_name, res1_id, res2_name, res2_id, distance)
    """
    residues = [r for r in chain.get_residues()
                if r.get_id()[0] == " " and "CA" in r][:max_residues]

    distances = []
    for i, r1 in enumerate(residues):
        for r2 in residues[i+1:]:
            d = ca_distance(r1, r2)
            distances.append((
                r1.get_resname(), r1.get_id()[1],
                r2.get_resname(), r2.get_id()[1],
                round(d, 3)
            ))
    return distances


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    if not os.path.exists(TIMP1_PDB):
        print(f"TIMP1 structure not found at {TIMP1_PDB}")
        print("Please run 03_fetch_pdb_structures.py first.")
        return

    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("TIMP1", TIMP1_PDB)

    # ── Active site analysis ──────────────────────────────────────────────────
    print("\n" + "=" * 55)
    print("TIMP1 ACTIVE SITE ANALYSIS (PDB: 1UEA)")
    print("=" * 55)
    print()
    print("TIMP1 inhibits MMPs by inserting its N-terminal Cys1 into the")
    print("MMP active site zinc. Let us see what surrounds Cys1:\n")

    nearby = find_nearby_residues(structure, chain_id="A", res_id=1, radius=5.0)

    print(f"Residues within 5 Å of Cys1 (chain A):")
    print(f"  {'Residue':<10} {'Number':<8} {'Notes'}")
    print(f"  {'-'*40}")

    active_site_data = []
    for res in nearby:
        resname = res.get_resname()
        resnum  = res.get_id()[1]
        note    = ""
        if resnum == 1:
            note = "← target residue (inhibitory tip)"
        active_site_data.append((resname, resnum, note))
        print(f"  {resname:<10} {resnum:<8} {note}")

    # ── Distance between Cys1 and Cys3 (disulfide bond pair) ─────────────────
    print()
    print("Key distances (Cys1 to neighbours):")
    chain_a  = structure[0]["A"]
    try:
        for partner_id in [3, 9, 13]:
            if partner_id in [r.get_id()[1] for r in chain_a.get_residues()]:
                d = ca_distance(chain_a[1], chain_a[partner_id])
                print(f"  Cys1 → residue {partner_id}: {d:.2f} Å")
    except KeyError:
        pass

    # ── Save active site to CSV ───────────────────────────────────────────────
    csv_path = os.path.join(RESULTS_DIR, "timp1_active_site.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["residue_name", "residue_number", "notes"])
        writer.writerows(active_site_data)
    print(f"\nActive site residues saved → {csv_path}")

    print("\n" + "=" * 55)
    print("Tip: Open pdb1uea.ent in PyMOL and run:")
    print('  select active_site, resi 1+3+9+13+68')
    print("  show sticks, active_site")
    print("  color yellow, active_site")


if __name__ == "__main__":
    main()
