"""
diagram_03_biology_concepts.py
================================
Generates 3 biological concept diagrams in plain English:
  A) What is liver fibrosis? (cascade diagram)
  B) How does AlphaFold work? (simplified pipeline)
  C) What does pLDDT mean? (colour-coded explanation)

Run: python diagrams/diagram_03_biology_concepts.py
Output: diagrams/concept_A_liver_fibrosis.png
        diagrams/concept_B_alphafold_explained.png
        diagrams/concept_C_plddt_explained.png
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

DARK_BLUE   = "#1A3A5C"
MID_BLUE    = "#2E6DA4"
LIGHT_BLUE  = "#AED6F1"
TEAL        = "#17A589"
LIGHT_TEAL  = "#A9DFBF"
ORANGE      = "#E67E22"
LIGHT_ORANGE= "#FAD7A0"
RED         = "#C0392B"
LIGHT_RED   = "#F1948A"
GREY        = "#ECF0F1"
DARK_GREY   = "#5D6D7E"
WHITE       = "#FFFFFF"
PURPLE      = "#7D3C98"
LIGHT_PURPLE= "#D7BDE2"
GREEN       = "#1E8449"
LIGHT_GREEN = "#A9DFBF"
GOLD        = "#D4AC0D"


# ══════════════════════════════════════════════════════════════════════════════
# DIAGRAM A — What is Liver Fibrosis?
# ══════════════════════════════════════════════════════════════════════════════
fig_a, ax_a = plt.subplots(figsize=(18, 13), facecolor=GREY)
ax_a.set_xlim(0, 18)
ax_a.set_ylim(0, 13)
ax_a.axis("off")

ax_a.text(9, 12.5, "What is Liver Fibrosis?",
          ha="center", fontsize=22, fontweight="bold", color=DARK_BLUE)
ax_a.text(9, 12.0,
          "A plain-English guide to the disease — and the proteins that drive it",
          ha="center", fontsize=12, color=DARK_GREY)

# --- Step boxes across the top ------------------------------------------------
cascade = [
    ("Liver Injury",      "Alcohol, viral\nhepatitis, fatty\nliver disease",   LIGHT_RED,    RED),
    ("Inflammation",      "Immune cells\nrush in. Alarm\nsignals released",    LIGHT_ORANGE, ORANGE),
    ("TGF-β1 Released",   "The master\nfibrosis trigger\nis turned on",        LIGHT_PURPLE, PURPLE),
    ("HSC Activation",    "Stellate cells\n(HSCs) wake up\n& multiply",        LIGHT_BLUE,   MID_BLUE),
    ("Collagen Made",     "COL1A1 builds\nscar matrix.\nTIMP1 stops\ncleanup", LIGHT_TEAL,   TEAL),
    ("Cirrhosis",         "Normal liver\ncells replaced\nby scar tissue",      "#F9E79F",    GOLD),
]

xs = np.linspace(1.3, 16.7, len(cascade))
for i, (title, body, fc, ec) in enumerate(cascade):
    box = FancyBboxPatch((xs[i]-1.1, 9.2), 2.2, 2.3,
                         boxstyle="round,pad=0.1",
                         facecolor=fc, edgecolor=ec, linewidth=2.5)
    ax_a.add_patch(box)
    ax_a.text(xs[i], 11.05, title, ha="center", va="center",
              fontsize=11, fontweight="bold", color=DARK_BLUE)
    ax_a.text(xs[i], 10.15, body, ha="center", va="center",
              fontsize=9.5, color=DARK_GREY, linespacing=1.4)
    if i < len(cascade) - 1:
        ax_a.annotate("", xy=(xs[i+1]-1.1-0.05, 10.35),
                      xytext=(xs[i]+1.1+0.05, 10.35),
                      arrowprops=dict(arrowstyle="-|>", color=DARK_BLUE,
                                      lw=2, mutation_scale=18))

# --- The four proteins explained ----------------------------------------------
ax_a.text(9, 8.7, "The Four Proteins We Study — and Why Each Matters",
          ha="center", fontsize=14, fontweight="bold", color=DARK_BLUE)

proteins = [
    ("TGF-β1", "P01137 | PDB: 1KTZ",
     "The ALARM button.\nWhen the liver is hurt,\nTGF-β1 is released\nand tells stellate\ncells to make scar.",
     LIGHT_RED, RED),
    ("ACTA2\n(α-SMA)", "P62736 | PDB: 3LUE",
     "The ACTIVATION marker.\nWhen stellate cells\nswitch on, they fill\nwith ACTA2 fibres and\nbecome contractile.",
     LIGHT_BLUE, MID_BLUE),
    ("COL1A1", "P02452 | PDB: 1BKV",
     "The SCAR BUILDER.\nCollagen Type I is the\nmain protein in scar\ntissue — tough, rigid,\nhard to break down.",
     LIGHT_TEAL, TEAL),
    ("TIMP1", "P01033 | PDB: 1UEA",
     "The ANTI-CLEANUP agent.\nBlocks the enzymes\n(MMPs) that would\nnormally dissolve\nextra scar tissue.",
     LIGHT_ORANGE, ORANGE),
]

pxs = [2.25, 6.25, 10.25, 14.25]
for i, (pname, pid, pdesc, fc, ec) in enumerate(proteins):
    box = FancyBboxPatch((pxs[i]-1.8, 1.2), 3.6, 6.8,
                         boxstyle="round,pad=0.12",
                         facecolor=fc, edgecolor=ec, linewidth=2.2)
    ax_a.add_patch(box)
    ax_a.text(pxs[i], 7.55, pname, ha="center", fontsize=13,
              fontweight="bold", color=DARK_BLUE, linespacing=1.2)
    ax_a.text(pxs[i], 6.95, pid, ha="center", fontsize=8.5,
              color=DARK_GREY, fontstyle="italic")
    ax_a.text(pxs[i], 5.2, pdesc, ha="center", va="center",
              fontsize=9.5, color=DARK_GREY, linespacing=1.5)
    # "drug target" badge
    badge_box = FancyBboxPatch((pxs[i]-1.3, 1.35), 2.6, 0.55,
                               boxstyle="round,pad=0.07",
                               facecolor=ec, edgecolor=ec, linewidth=0)
    ax_a.add_patch(badge_box)
    ax_a.text(pxs[i], 1.63, "Drug Target", ha="center",
              fontsize=9, fontweight="bold", color=WHITE)

plt.tight_layout()
plt.savefig("diagrams/concept_A_liver_fibrosis.png", dpi=150,
            bbox_inches="tight", facecolor=GREY)
plt.close()
print("✅  Saved: diagrams/concept_A_liver_fibrosis.png")


# ══════════════════════════════════════════════════════════════════════════════
# DIAGRAM B — How AlphaFold Works
# ══════════════════════════════════════════════════════════════════════════════
fig_b, ax_b = plt.subplots(figsize=(18, 12), facecolor=GREY)
ax_b.set_xlim(0, 18)
ax_b.set_ylim(0, 12)
ax_b.axis("off")

ax_b.text(9, 11.5, "How Does AlphaFold Work?",
          ha="center", fontsize=22, fontweight="bold", color=DARK_BLUE)
ax_b.text(9, 11.0,
          "A simple explanation of the AI system that predicts protein shapes from sequence alone",
          ha="center", fontsize=12, color=DARK_GREY)

# Pipeline steps
af_steps = [
    ("INPUT",
     "Amino Acid Sequence",
     "You give AlphaFold a string\nof letters — the protein sequence.\nE.g. MHIKL... (one letter per\namino acid)",
     LIGHT_BLUE, MID_BLUE, "1"),
    ("SEARCH",
     "Find Related Sequences\n(Multiple Sequence Alignment)",
     "AlphaFold searches millions of\nsequences from other species.\nAmino acids that always change\ntogether → probably touching in 3D!",
     LIGHT_PURPLE, PURPLE, "2"),
    ("AI MODEL",
     "Evoformer Neural Network\n(48 layers of attention)",
     "The deep learning model processes\nboth sequence + co-evolution patterns.\nLike reading between the lines\nof evolution to infer 3D structure.",
     LIGHT_ORANGE, ORANGE, "3"),
    ("OUTPUT",
     "3D Coordinates + pLDDT Scores",
     "AlphaFold outputs x,y,z positions\nfor every atom in the protein,\nplus a confidence score (pLDDT)\nfor each residue.",
     LIGHT_GREEN, GREEN, "4"),
]

step_xs = [2.1, 6.3, 11.3, 15.9]
step_ws = [3.8, 4.2, 4.5, 3.8]

for i, (tag, title, body, fc, ec, num) in enumerate(af_steps):
    w = step_ws[i]
    box = FancyBboxPatch((step_xs[i]-w/2, 4.5), w, 5.0,
                         boxstyle="round,pad=0.12",
                         facecolor=fc, edgecolor=ec, linewidth=2.5)
    ax_b.add_patch(box)

    # Number badge
    badge = plt.Circle((step_xs[i]-w/2+0.3, 9.2), 0.25, color=ec, zorder=5)
    ax_b.add_patch(badge)
    ax_b.text(step_xs[i]-w/2+0.3, 9.2, num,
              ha="center", va="center", fontsize=11, fontweight="bold",
              color=WHITE, zorder=6)

    ax_b.text(step_xs[i], 9.15, tag,
              ha="center", fontsize=9, fontstyle="italic",
              color=DARK_GREY, fontfamily="monospace")
    ax_b.text(step_xs[i], 8.5, title,
              ha="center", va="center", fontsize=11, fontweight="bold",
              color=DARK_BLUE, linespacing=1.3)
    ax_b.text(step_xs[i], 6.5, body,
              ha="center", va="center", fontsize=9.5,
              color=DARK_GREY, linespacing=1.5)

    # Arrow between steps
    if i < len(af_steps) - 1:
        next_x = step_xs[i+1]
        next_w = step_ws[i+1]
        ax_b.annotate("", xy=(next_x - next_w/2 - 0.05, 7.0),
                      xytext=(step_xs[i] + w/2 + 0.05, 7.0),
                      arrowprops=dict(arrowstyle="-|>", color=DARK_BLUE,
                                      lw=2.5, mutation_scale=20))

# Key insight box
insight = FancyBboxPatch((1.0, 1.0), 16, 3.0,
                          boxstyle="round,pad=0.15",
                          facecolor="#EAF2FF", edgecolor=MID_BLUE, linewidth=2)
ax_b.add_patch(insight)
ax_b.text(9, 3.6, "Key Insight — Why AlphaFold is Revolutionary",
          ha="center", fontsize=12, fontweight="bold", color=MID_BLUE)
ax_b.text(9, 2.9,
          "Before AlphaFold: determining a protein's 3D structure required months in a lab,\n"
          "expensive equipment (X-ray crystallography or cryo-EM), and often failed completely.",
          ha="center", va="center", fontsize=10, color=DARK_GREY, linespacing=1.5)
ax_b.text(9, 1.8,
          "After AlphaFold: anyone with a sequence and a laptop can get a high-quality\n"
          "structure prediction in minutes — for FREE. It opened structural biology to the whole world.",
          ha="center", va="center", fontsize=10, color=DARK_GREY, linespacing=1.5)

plt.tight_layout()
plt.savefig("diagrams/concept_B_alphafold_explained.png", dpi=150,
            bbox_inches="tight", facecolor=GREY)
plt.close()
print("✅  Saved: diagrams/concept_B_alphafold_explained.png")


# ══════════════════════════════════════════════════════════════════════════════
# DIAGRAM C — What does pLDDT mean?
# ══════════════════════════════════════════════════════════════════════════════
fig_c, ax_c = plt.subplots(figsize=(18, 13), facecolor=GREY)
ax_c.set_xlim(0, 18)
ax_c.set_ylim(0, 13)
ax_c.axis("off")

ax_c.text(9, 12.5, "What Does pLDDT Mean?",
          ha="center", fontsize=22, fontweight="bold", color=DARK_BLUE)
ax_c.text(9, 12.0,
          "AlphaFold's confidence score — and why low confidence is not a failure",
          ha="center", fontsize=12, color=DARK_GREY)

# Colour scale bar
bar_colors = ["#FF7D45", "#FFDB13", "#65CBF3", "#006EFF"]
bar_labels  = ["0", "50", "70", "90", "100"]
bar_xs = [1.5, 5.5, 10.0, 14.0]
bar_w  = 4.0
bar_y  = 9.0

for i, (fc) in enumerate(bar_colors):
    box = FancyBboxPatch((bar_xs[i], bar_y), bar_w, 1.6,
                         boxstyle="square,pad=0",
                         facecolor=fc, edgecolor=WHITE, linewidth=1.5)
    ax_c.add_patch(box)
ax_c.text(1.5,  bar_y - 0.35, "0",   ha="center", fontsize=10, color=DARK_GREY)
ax_c.text(5.5,  bar_y - 0.35, "50",  ha="center", fontsize=10, color=DARK_GREY)
ax_c.text(10.0, bar_y - 0.35, "70",  ha="center", fontsize=10, color=DARK_GREY)
ax_c.text(14.0, bar_y - 0.35, "90",  ha="center", fontsize=10, color=DARK_GREY)
ax_c.text(18.0, bar_y - 0.35, "100", ha="center", fontsize=10, color=DARK_GREY)
ax_c.text(9, bar_y + 2.0, "pLDDT Score  (0 = no confidence   →   100 = perfect confidence)",
          ha="center", fontsize=11, color=DARK_GREY)

# Four explanation cards
cards = [
    ("< 50",
     "Very Low Confidence\n(Orange / Red)",
     "This region is probably\nINTRINSICALLY DISORDERED\n— it doesn't have a\nfixed shape in real life.\n\nNot a prediction failure!\nDisorder is biologically real.",
     "#FF7D45", WHITE),
    ("50 – 70",
     "Low Confidence\n(Yellow)",
     "Flexible region.\nThe overall shape is\nprobably roughly right,\nbut don't trust the\nexact atom positions.\n\nOften: surface loops,\nlinkers between domains.",
     "#FFDB13", DARK_BLUE),
    ("70 – 90",
     "Good Confidence\n(Light Blue)",
     "Well-structured region.\nThe backbone shape is\nlikely correct.\nSide-chain positions\nmay have some uncertainty.\n\nGood for drug design!",
     "#65CBF3", DARK_BLUE),
    ("> 90",
     "Very High Confidence\n(Dark Blue)",
     "Rigid, well-folded domain.\nAtom positions are as\naccurate as an X-ray\ncrystal structure.\n\nTrust every detail here.",
     "#006EFF", WHITE),
]

card_xs = [2.25, 6.25, 10.5, 14.75]
for i, (score, title, body, fc, tc) in enumerate(cards):
    box = FancyBboxPatch((card_xs[i]-1.85, 1.2), 3.7, 7.1,
                         boxstyle="round,pad=0.12",
                         facecolor=fc, edgecolor=DARK_BLUE, linewidth=2.2)
    ax_c.add_patch(box)
    ax_c.text(card_xs[i], 7.95, score,
              ha="center", fontsize=18, fontweight="bold", color=tc)
    ax_c.text(card_xs[i], 7.1, title,
              ha="center", va="center", fontsize=11, fontweight="bold",
              color=tc, linespacing=1.3)
    ax_c.text(card_xs[i], 4.5, body,
              ha="center", va="center", fontsize=9.5,
              color=tc if tc == WHITE else DARK_GREY, linespacing=1.5)

# Footer insight
insight = FancyBboxPatch((1.0, 0.15), 16, 0.85,
                          boxstyle="round,pad=0.1",
                          facecolor=LIGHT_TEAL, edgecolor=TEAL, linewidth=2)
ax_c.add_patch(insight)
ax_c.text(9, 0.57,
          "Key insight: In Project 2, LOXL2's propeptide region (residues 315–489) "
          "scores below 50 — and that is correct! It is disordered by design.",
          ha="center", va="center", fontsize=10, color=DARK_BLUE)

plt.tight_layout()
plt.savefig("diagrams/concept_C_plddt_explained.png", dpi=150,
            bbox_inches="tight", facecolor=GREY)
plt.close()
print("✅  Saved: diagrams/concept_C_plddt_explained.png")
