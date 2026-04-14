"""
diagram_01_project1_workflow.py
================================
Generates a step-by-step workflow flowchart for Project 1.
Shows what each script does in plain English, with arrows connecting each step.

Run: python diagrams/diagram_01_project1_workflow.py
Output: diagrams/project1_workflow.png
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ── Colour palette ─────────────────────────────────────────────────────────────
DARK_BLUE   = "#1A3A5C"
MID_BLUE    = "#2E6DA4"
LIGHT_BLUE  = "#AED6F1"
TEAL        = "#17A589"
LIGHT_TEAL  = "#A9DFBF"
ORANGE      = "#E67E22"
LIGHT_ORANGE= "#FAD7A0"
GREY        = "#ECF0F1"
DARK_GREY   = "#5D6D7E"
WHITE       = "#FFFFFF"

fig = plt.figure(figsize=(18, 22), facecolor=GREY)
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 18)
ax.set_ylim(0, 22)
ax.axis("off")

# ── Title ──────────────────────────────────────────────────────────────────────
ax.text(9, 21.3, "Project 1 — Workflow Flowchart",
        ha="center", va="center", fontsize=22, fontweight="bold",
        color=DARK_BLUE, fontfamily="DejaVu Sans")
ax.text(9, 20.8, "Protein Structure Visualisation of Liver Cirrhosis-Linked Proteins",
        ha="center", va="center", fontsize=13, color=DARK_GREY)

# ── Helper: draw a step box ─────────────────────────────────────────────────────
def draw_step(ax, x, y, w, h, step_num, script_name, title, bullets,
              box_color, text_color=DARK_BLUE, icon=""):
    """Draws a rounded rectangle representing one pipeline step."""
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.15",
                         facecolor=box_color, edgecolor=DARK_BLUE, linewidth=2)
    ax.add_patch(box)

    # Step badge
    badge = plt.Circle((x - w/2 + 0.35, y + h/2 - 0.25), 0.22,
                        color=DARK_BLUE, zorder=5)
    ax.add_patch(badge)
    ax.text(x - w/2 + 0.35, y + h/2 - 0.25, str(step_num),
            ha="center", va="center", fontsize=10, fontweight="bold",
            color=WHITE, zorder=6)

    # Script filename
    ax.text(x, y + h/2 - 0.22, script_name,
            ha="center", va="center", fontsize=9, color=DARK_GREY,
            fontstyle="italic", fontfamily="monospace")

    # Title
    ax.text(x, y + h/2 - 0.55, title,
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=text_color)

    # Bullet points
    for i, b in enumerate(bullets):
        ax.text(x - w/2 + 0.3, y + h/2 - 0.90 - i * 0.38, f"• {b}",
                ha="left", va="center", fontsize=9.5, color=DARK_GREY)


# ── Helper: draw arrow ─────────────────────────────────────────────────────────
def draw_arrow(ax, x, y_start, y_end, label=""):
    ax.annotate("", xy=(x, y_end + 0.12), xytext=(x, y_start - 0.12),
                arrowprops=dict(arrowstyle="-|>", color=MID_BLUE,
                                lw=2.5, mutation_scale=20))
    if label:
        ax.text(x + 0.25, (y_start + y_end) / 2, label,
                ha="left", va="center", fontsize=9, color=DARK_GREY, fontstyle="italic")


# ── Step definitions ───────────────────────────────────────────────────────────
steps = [
    {
        "y": 19.1, "h": 1.5,
        "step_num": 1,
        "script_name": "01_fetch_sequences.py",
        "title": "Download Protein Sequences from UniProt",
        "bullets": [
            "Connects to UniProt — the world's protein database",
            "Downloads sequences for TGF-β1, ACTA2, COL1A1, TIMP1",
            "Saves each sequence as a .fasta file (standard text format)",
        ],
        "box_color": LIGHT_BLUE,
        "arrow_label": "4 FASTA files ready",
    },
    {
        "y": 16.85, "h": 1.5,
        "step_num": 2,
        "script_name": "02_protparam_analysis.py",
        "title": "Calculate Physico-Chemical Properties",
        "bullets": [
            "Reads each FASTA file and runs ProteinAnalysis",
            "Calculates: molecular weight, charge (pI), hydrophobicity (GRAVY)",
            "Estimates secondary structure: % helix, % sheet, % coil",
        ],
        "box_color": LIGHT_TEAL,
        "arrow_label": "property table + charts",
    },
    {
        "y": 14.6, "h": 1.5,
        "step_num": 3,
        "script_name": "03_fetch_pdb_structures.py",
        "title": "Download 3D Structures from the PDB",
        "bullets": [
            "Fetches experimental 3D structure files (.pdb) from RCSB",
            "Parses the structure hierarchy: protein → chain → residue → atom",
            "Reports how many residues and atoms are in each structure",
        ],
        "box_color": LIGHT_BLUE,
        "arrow_label": "4 .pdb files parsed",
    },
    {
        "y": 12.35, "h": 1.5,
        "step_num": 4,
        "script_name": "04_residue_distances.py",
        "title": "Measure Distances Between Key Residues",
        "bullets": [
            "Uses NeighborSearch to find residues within 5 Ångstroms of each other",
            "Calculates Cα–Cα distances between specific residue pairs",
            "Identifies active site neighbourhood (e.g. TIMP1 Cys1 inhibitory tip)",
        ],
        "box_color": LIGHT_TEAL,
        "arrow_label": "distance tables saved",
    },
    {
        "y": 10.1, "h": 1.5,
        "step_num": 5,
        "script_name": "PyMOL: visualise_tgfb1.pml",
        "title": "Visualise Individual Protein in PyMOL",
        "bullets": [
            "Loads the PDB file into PyMOL (the molecular 'microscope')",
            "Colours by secondary structure: red=helix, yellow=sheet, green=coil",
            "Highlights the receptor-binding loop; adds semi-transparent surface",
        ],
        "box_color": LIGHT_ORANGE,
        "arrow_label": "PNG image rendered",
    },
    {
        "y": 7.85, "h": 1.5,
        "step_num": 6,
        "script_name": "PyMOL: compare_all_proteins.pml",
        "title": "Compare All Four Proteins Side-by-Side",
        "bullets": [
            "Loads TGF-β1, ACTA2, COL1A1, and TIMP1 simultaneously",
            "Gives each a distinct colour; arranges them side-by-side",
            "Saves the session (.pse) and renders a comparison image",
        ],
        "box_color": LIGHT_ORANGE,
        "arrow_label": "comparison figure saved",
    },
]

# ── Draw steps ─────────────────────────────────────────────────────────────────
BOX_W = 13
for i, s in enumerate(steps):
    draw_step(ax, 9, s["y"], BOX_W, s["h"],
              s["step_num"], s["script_name"], s["title"],
              s["bullets"], s["box_color"])
    if i < len(steps) - 1:
        draw_arrow(ax, 9,
                   s["y"] - s["h"]/2,
                   steps[i+1]["y"] + steps[i+1]["h"]/2,
                   s.get("arrow_label", ""))

# ── Output box at the bottom ───────────────────────────────────────────────────
out_box = FancyBboxPatch((2.5, 5.5), 13, 1.6,
                         boxstyle="round,pad=0.15",
                         facecolor="#D5F5E3", edgecolor=TEAL, linewidth=2.5)
ax.add_patch(out_box)
ax.text(9, 6.6, "Final Outputs",
        ha="center", va="center", fontsize=13, fontweight="bold", color=TEAL)
ax.text(9, 6.1, "ProtParam summary table  •  Secondary structure bar charts  •  "
        "3D protein images (PNG)  •  Comparison figure  •  Saved PyMOL session",
        ha="center", va="center", fontsize=10, color=DARK_GREY)

# ── Legend: what the colours mean ─────────────────────────────────────────────
legend_items = [
    mpatches.Patch(color=LIGHT_BLUE,   label="BioPython script (Python)"),
    mpatches.Patch(color=LIGHT_TEAL,   label="BioPython analysis step"),
    mpatches.Patch(color=LIGHT_ORANGE, label="PyMOL visualisation script"),
    mpatches.Patch(color="#D5F5E3",    label="Final output / result"),
]
ax.legend(handles=legend_items, loc="lower center",
          bbox_to_anchor=(0.5, 0.01), ncol=4,
          fontsize=10, framealpha=0.9,
          edgecolor=DARK_BLUE, facecolor=WHITE)

plt.savefig("diagrams/project1_workflow.png", dpi=150,
            bbox_inches="tight", facecolor=GREY)
plt.close()
print("✅  Saved: diagrams/project1_workflow.png")
