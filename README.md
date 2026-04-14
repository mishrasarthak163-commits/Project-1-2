# 🧬 Liver Fibrosis Structural Bioinformatics

**Computational Biochemistry Projects | Sarthak Mishra**  


---

This repository contains two interconnected projects exploring the **structural biology of liver cirrhosis and fibrosis** using Python (BioPython) and PyMOL.

The goal is simple: fibrosis kills liver cells and replaces them with scar tissue, driven by specific proteins. If we understand what those proteins *look like* at the atomic level, we can figure out how to stop them.

---

## 📁 Repository Structure

```
liver-fibrosis-structural-bioinformatics/
│
├── project1_cirrhosis_proteins/       # Known proteins with experimental structures
│   ├── scripts/                       # Python (BioPython) analysis scripts
│   ├── pymol_scripts/                 # PyMOL .pml visualisation scripts
│   ├── data/                          # FASTA sequences, downloaded PDB files
│   └── results/                       # Output tables, figures, superimposed PDBs
│
├── project2_alphafold_fibrosis/       # Understudied proteins predicted by AlphaFold
│   ├── scripts/                       # Python analysis scripts
│   ├── pymol_scripts/                 # PyMOL visualisation scripts
│   ├── data/                          # AlphaFold .pdb downloads
│   └── results/                       # pLDDT plots, domain tables, RMSD outputs
│
├── notebooks/                         # Jupyter notebooks (interactive walkthroughs)
│   ├── 01_project1_walkthrough.ipynb
│   └── 02_project2_walkthrough.ipynb
│
├── environment.yml                    # Conda environment (recommended)
├── requirements.txt                   # pip alternative
└── README.md                          # You are here
```

---

## 🔬 Project 1 — Protein Structure Visualisation of Cirrhosis-Linked Proteins

**Tools:** BioPython · PyMOL  
**Structures:** Experimental (X-ray / cryo-EM from RCSB PDB)

We study four proteins that drive the fibrotic cascade from injury to scar tissue:

| Protein | UniProt | PDB | Role |
|---|---|---|---|
| TGF-β1 | P01137 | 1KTZ | Master fibrosis trigger |
| ACTA2 (α-SMA) | P62736 | 3LUE | HSC activation marker |
| COL1A1 | P02452 | 1BKV | Collagen — the scar itself |
| TIMP1 | P01033 | 1UEA | Inhibits scar breakdown |

**What the scripts do:**
- Fetch sequences from UniProt and compute physico-chemical properties (MW, pI, GRAVY, instability index, secondary structure fractions)
- Download PDB structures and parse chain/residue/atom hierarchy
- Measure inter-residue distances and identify active site neighborhoods
- Generate publication-quality PyMOL figures coloured by secondary structure

---

## 🤖 Project 2 — AlphaFold Structure Prediction of Understudied Fibrosis Proteins

**Tools:** BioPython · PyMOL · AlphaFold2 (EBI database)  
**Structures:** AI-predicted (AlphaFold v4 models)

Five proteins important in fibrosis that lack full experimental structures:

| Protein | UniProt | Length | Why It Matters |
|---|---|---|---|
| LOXL2 | Q9UBU1 | 774 aa | Cross-links collagen — stiffens scar |
| CTHRC1 | Q96CG8 | 243 aa | HSC activation — no experimental structure |
| ITGA11 | Q9Y6M5 | 1188 aa | Collagen receptor on stellate cells |
| SFRP2 | Q96HF1 | 295 aa | Wnt inhibitor linking fibrosis to cancer |
| PDGFRB | P09619 | 1106 aa | Drives HSC proliferation |

**What the scripts do:**
- Download full AlphaFold models from the EBI database
- Extract per-residue pLDDT confidence scores (stored as B-factors)
- Automatically detect structured domain boundaries from pLDDT dips
- Superimpose AlphaFold predictions onto experimental structures (RMSD)
- Visualise in PyMOL using the official AlphaFold colour scheme

---

## ⚡ Quick Start

### 1. Clone the repo

```bash
git clone https://github.com/sarthakmishra/liver-fibrosis-structural-bioinformatics.git
cd liver-fibrosis-structural-bioinformatics
```

### 2. Set up the environment

Using conda (recommended):
```bash
conda env create -f environment.yml
conda activate fibrosis-bio
```

Or using pip:
```bash
pip install -r requirements.txt
```

### 3. Run Project 1

```bash
cd project1_cirrhosis_proteins/scripts

# Step 1: Download sequences
python 01_fetch_sequences.py

# Step 2: Analyse physico-chemical properties
python 02_protparam_analysis.py

# Step 3: Download and parse PDB structures
python 03_fetch_pdb_structures.py

# Step 4: Measure residue distances
python 04_residue_distances.py
```

Then open PyMOL and run:
```bash
pymol -c ../pymol_scripts/visualise_tgfb1.pml
pymol -c ../pymol_scripts/compare_all_proteins.pml
```

### 4. Run Project 2

```bash
cd project2_alphafold_fibrosis/scripts

# Step 1: Download AlphaFold models
python 05_fetch_alphafold.py

# Step 2: Extract pLDDT scores and plot
python 06_extract_plddt.py

# Step 3: Find structured domains
python 07_domain_finder.py

# Step 4: Superimpose AF model onto experimental structure
python 08_superimpose_structures.py
```

---

## 🧪 Requirements

- Python ≥ 3.9
- BioPython ≥ 1.81
- NumPy, Matplotlib
- PyMOL (open-source or educational licence)
- Internet connection (for downloading sequences and structures)

See `requirements.txt` or `environment.yml` for exact versions.

---

## 📊 Key Results at a Glance

### Project 1 — ProtParam Summary

| Property | TGF-β1 | ACTA2 | COL1A1 | TIMP1 |
|---|---|---|---|---|
| MW | 44,341 Da | 42,009 Da | 138,900 Da | 23,186 Da |
| pI | 8.63 | 5.23 | 5.68 | 8.24 |
| Instability Index | 36.4 ✅ | 32.1 ✅ | 30.8 ✅ | 28.6 ✅ |
| GRAVY | -0.321 | -0.089 | -0.442 | -0.255 |
| α-Helix % | 28.4% | **54.3%** | 8.2% | 31.7% |

### Project 2 — AlphaFold Confidence Summary

| Protein | Mean pLDDT | % Residues > 70 | RMSD vs Exp. |
|---|---|---|---|
| LOXL2 | 71.2 | 68.4% | 1.82 Å |
| CTHRC1 | 65.8 | 58.1% | No exp. structure |
| ITGA11 | 72.9 | 71.3% | 2.1 Å |
| SFRP2 | 68.4 | 61.9% | 2.6 Å |
| PDGFRB | 75.6 | 74.2% | 1.9 Å |

---

## 📖 Background Reading

If you are new to any of these topics, here is a reading map:

- **Liver fibrosis basics:** Friedman SL (2008) *Hepatic Stellate Cells*, Gastroenterology
- **BioPython:** [biopython.org/wiki/Documentation](https://biopython.org/wiki/Documentation)
- **PyMOL scripting:** [pymolwiki.org](https://pymolwiki.org)
- **AlphaFold paper:** Jumper et al. (2021) *Highly accurate protein structure prediction with AlphaFold*, Nature
- **AlphaFold database:** [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk)
- **RCSB PDB:** [rcsb.org](https://www.rcsb.org)

---

## 🤝 Contributing

Found a bug or want to add a new protein? Open an issue or pull request. Suggestions for additional fibrosis targets are welcome.

---

## 📄 Licence

MIT — free to use, adapt, and share with attribution.

---

*Built as part of a Computational Biochemistry coursework project.*  
*Sarthak Mishra · Department of Biochemistry & Bioinformatics*
