"""
01_fetch_sequences.py
─────────────────────
Downloads amino acid sequences for the four cirrhosis-linked proteins
from the UniProt database and saves them as FASTA files.

Proteins:
  TGF-β1   (P01137) – master fibrosis trigger
  ACTA2    (P62736) – HSC activation marker
  COL1A1   (P02452) – main collagen in scar tissue
  TIMP1    (P01033) – inhibitor of scar breakdown

Usage:
  python 01_fetch_sequences.py

Output:
  ../data/TGF_beta1.fasta
  ../data/ACTA2.fasta
  ../data/COL1A1.fasta
  ../data/TIMP1.fasta
"""

import urllib.request
import time
import os

# ── Config ────────────────────────────────────────────────────────────────────

OUTPUT_DIR = "../data"

# UniProt accession IDs for our four proteins
PROTEINS = {
    "TGF_beta1": "P01137",
    "ACTA2":     "P62736",
    "COL1A1":    "P02452",
    "TIMP1":     "P01033",
}

# ── Main ──────────────────────────────────────────────────────────────────────

def fetch_uniprot_fasta(accession: str) -> str:
    """
    Downloads a protein's FASTA sequence from UniProt.

    UniProt gives each protein a unique accession number — think of it
    like a barcode. We build a URL with that barcode and download the
    sequence as plain text.

    Args:
        accession: UniProt accession ID (e.g. 'P01137')

    Returns:
        FASTA-formatted string
    """
    url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
    with urllib.request.urlopen(url) as response:
        return response.read().decode("utf-8")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Downloading protein sequences from UniProt...\n")

    for name, accession in PROTEINS.items():
        print(f"  [{accession}] {name} ...", end=" ", flush=True)

        fasta = fetch_uniprot_fasta(accession)

        out_path = os.path.join(OUTPUT_DIR, f"{name}.fasta")
        with open(out_path, "w") as f:
            f.write(fasta)

        # Count residues (FASTA lines starting with '>' are headers, skip them)
        seq_lines = [l for l in fasta.splitlines() if not l.startswith(">")]
        n_residues = sum(len(l) for l in seq_lines)

        print(f"done  ({n_residues} residues → {out_path})")

        # Be polite to the server — one request per second
        time.sleep(1)

    print("\nAll sequences downloaded successfully!")
    print(f"Files saved to: {os.path.abspath(OUTPUT_DIR)}")


if __name__ == "__main__":
    main()
