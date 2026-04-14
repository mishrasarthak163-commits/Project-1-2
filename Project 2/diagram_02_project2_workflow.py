"""
diagram_02_project2_workflow.py
================================
Generates a step-by-step workflow flowchart for Project 2.
Shows what each AlphaFold script does in plain English.

Run: python diagrams/diagram_02_project2_workflow.py
Output: diagrams/project2_workflow.png
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

DARK_BLUE    = "#1A3A5C"
MID_BLUE     = "#2E6DA4"
LIGHT_BLUE   = "#AED6F1"
PURPLE       = "#7D3C98"
LIGHT_PURPLE = "#D7BDE2"
ORANGE       = "#E67E22"
LIGHT_ORANGE = "#FAD7A0"
GREEN        = "#1E8449"
LIGHT_GREEN  = "#A9DFBF"
GREY         = "#ECF0F1"
DARK_GREY    = "#5D6D7E"
WHITE        = "#FFFFFF"
GOLD         = "#F1C40F"

fig = plt.figure(figsize=(18, 24), facecolor=GREY)
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 18)
ax.set_ylim(0, 24)
ax.axis("off")

# Title
ax.text(9, 23.3, "Project 2 — Workflow Flowchart",
        ha="center", va="center", fontsize=22, fontweight="bold",
        color=DARK_BLUE)
ax.text(9, 22.8, "AlphaFold Structure Prediction of Understudied Fibrosis Proteins",
        ha="center", va="center", fontsize=13, color=DARK_GREY)

def draw_step(ax, x, y, w, h, step_num, script_name, title, bullets,
              box_color, text_color=DARK_BLUE):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.15",
                         facecolor=box_color, edgecolor=DARK_BLUE, linewidth=2)
    ax.add_patch(box)
    badge = plt.Circle((x - w/2 + 0.35, y + h/2 - 0.25), 0.22,
                        color=DARK_BLUE, zorder=5)
    ax.add_patch(badge)
    ax.text(x - w/2 + 0.35, y + h/2 - 0.25, str(step_num),
            ha="center", va="center", fontsize=10, fontweight="bold",
            color=WHITE, zorder=6)
    ax.text(x, y + h/2 - 0.22, script_name,
            ha="center", va="center", fontsize=9, color=DARK_GREY,
            fontstyle="italic", fontfamily="monospace")
    ax.text(x, y + h/2 - 0.55, title,
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=text_color)
    for i, b in enumerate(bullets):
        ax.text(x - w/2 + 0.3, y + h/2 - 0.90 - i * 0.38, f"• {b}",
                ha="left", va="center", fontsize=9.5, color=DARK_GREY)

def draw_arrow(ax, x, y_start, y_end, label=""):
    ax.annotate("", xy=(x, y_end + 0.12), xytext=(x, y_start - 0.12),
                arrowprops=dict(arrowstyle="-|>", color=MID_BLUE,
                                lw=2.5, mutation_scale=20))
    if label:
        ax.text(x + 0.25, (y_start + y_end) / 2, label,
                ha="left", va="center", fontsize=9, color=DARK_GREY, fontstyle="italic")

steps = [
    {
        "y": 21.2, "h": 1.5,
        "step_num": 1,
        "script_name": "05_fetch_alphafold.py",
        "title": "Download AlphaFold Predicted Structures",
        "bullets": [
            "Connects to the AlphaFold EBI database (free, open access)",
            "Downloads full-length predicted .pdb files for 5 proteins",
            "Proteins: LOXL2, CTHRC1, ITGA11, SFRP2, PDGFRB",
        ],
        "box_color": LIGHT_PURPLE,
        "arrow_label": "5 AlphaFold .pdb files",
    },
    {
        "y": 18.95, "h": 1.5,
        "step_num": 2,
        "script_name": "06_extract_plddt.py",
        "title": "Extract pLDDT Confidence Scores",
        "bullets": [
            "Reads each AlphaFold PDB file — pLDDT is stored as the B-factor",
            "Extracts one confidence score (0–100) per residue",
            "Plots score along the sequence: where is AlphaFold confident?",
        ],
        "box_color": LIGHT_BLUE,
        "arrow_label": "pLDDT arrays + line plots",
    },
    {
        "y": 16.7, "h": 1.5,
        "step_num": 3,
        "script_name": "07_domain_finder.py",
        "title": "Automatically Detect Structured Domains",
        "bullets": [
            "Scans the pLDDT profile for regions consistently above 70",
            "Labels these as structured domains (≥ 30 residues long)",
            "Identifies flexible linkers where pLDDT dips below 50",
        ],
        "box_color": LIGHT_GREEN,
        "arrow_label": "domain boundary table",
    },
    {
        "y": 14.45, "h": 1.5,
        "step_num": 4,
        "script_name": "08_superimpose_structures.py",
        "title": "Compare AlphaFold vs Experimental Structure",
        "bullets": [
            "Loads both the AlphaFold model and the experimental PDB structure",
            "Overlays them in 3D space using Cα atom alignment",
            "Calculates RMSD — how far apart (in Ångstroms) the two models are",
        ],
        "box_color": LIGHT_ORANGE,
        "arrow_label": "RMSD value + superimposed .pdb",
    },
    {
        "y": 12.2, "h": 1.5,
        "step_num": 5,
        "script_name": "PyMOL: colour_by_plddt.pml",
        "title": "Visualise Structure Coloured by Confidence",
        "bullets": [
            "Opens AlphaFold model in PyMOL — the molecular visualisation tool",
            "Colours each residue by its pLDDT: blue=confident, orange=disordered",
            "Labels key domains; adds semi-transparent surface for context",
        ],
        "box_color": LIGHT_PURPLE,
        "arrow_label": "coloured PNG rendered",
    },
    {
        "y": 9.95, "h": 1.5,
        "step_num": 6,
        "script_name": "PyMOL: compare_experimental_vs_alphafold.pml",
        "title": "Side-by-Side Comparison Figure",
        "bullets": [
            "Loads both experimental (blue) and AlphaFold (orange) structures",
            "Places them side-by-side after superimposition",
            "Highlights extra regions covered by AlphaFold but not by experiment",
        ],
        "box_color": LIGHT_ORANGE,
        "arrow_label": "comparison figure saved",
    },
    {
        "y": 7.7, "h": 1.5,
        "step_num": 7,
        "script_name": "PyMOL: visualise_all_alphafold.pml",
        "title": "All Five Proteins Together",
        "bullets": [
            "Loads all five AlphaFold models simultaneously",
            "Arranges them side-by-side with distinct colours",
            "Saves the PyMOL session for interactive exploration",
        ],
        "box_color": LIGHT_PURPLE,
        "arrow_label": "5-protein panel figure",
    },
]

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

# Output box
out_box = FancyBboxPatch((2.5, 5.5), 13, 1.7,
                         boxstyle="round,pad=0.15",
                         facecolor="#D5F5E3", edgecolor=GREEN, linewidth=2.5)
ax.add_patch(out_box)
ax.text(9, 6.65, "Final Outputs",
        ha="center", va="center", fontsize=13, fontweight="bold", color=GREEN)
ax.text(9, 6.15, "pLDDT line plots  •  Domain boundary tables  •  RMSD values  •  "
        "pLDDT-coloured 3D images  •  Experimental vs AlphaFold comparison figures",
        ha="center", va="center", fontsize=10, color=DARK_GREY)

# Legend
legend_items = [
    mpatches.Patch(color=LIGHT_PURPLE, label="AlphaFold / prediction step"),
    mpatches.Patch(color=LIGHT_BLUE,   label="BioPython analysis step"),
    mpatches.Patch(color=LIGHT_GREEN,  label="Domain detection step"),
    mpatches.Patch(color=LIGHT_ORANGE, label="PyMOL visualisation / comparison"),
    mpatches.Patch(color="#D5F5E3",    label="Final output / result"),
]
ax.legend(handles=legend_items, loc="lower center",
          bbox_to_anchor=(0.5, 0.005), ncol=3,
          fontsize=10, framealpha=0.9,
          edgecolor=DARK_BLUE, facecolor=WHITE)

plt.savefig("diagrams/project2_workflow.png", dpi=150,
            bbox_inches="tight", facecolor=GREY)
plt.close()
print("✅  Saved: diagrams/project2_workflow.png")
