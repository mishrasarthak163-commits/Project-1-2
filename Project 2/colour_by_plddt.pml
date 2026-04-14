# colour_by_plddt.pml
# ─────────────────────────────────────────────────────────────────────────────
# Visualises the AlphaFold LOXL2 model (Q9UBU1) coloured by pLDDT confidence.
#
# Why colour by pLDDT?
#   In an AlphaFold PDB file, the B-factor column stores the pLDDT score
#   (0–100) instead of atomic vibration. Colouring by B-factor therefore
#   gives an instant visual map of confidence:
#
#   Dark blue  (>90) → very confident, well-structured
#   Light blue (70–90) → confident
#   Yellow     (50–70) → low confidence, flexible
#   Orange     (<50)  → likely disordered — no stable 3D structure
#
# How to run:
#   Place LOXL2_AF.pdb in the same directory as this script, then:
#   pymol -c colour_by_plddt.pml
#
# Output: LOXL2_pLDDT_coloured.png
# ─────────────────────────────────────────────────────────────────────────────

# ── Load structure ────────────────────────────────────────────────────────────
load LOXL2_AF.pdb, LOXL2

# ── Clean display ─────────────────────────────────────────────────────────────
hide everything
bg_color white
show cartoon, LOXL2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# ── Colour by B-factor (= pLDDT) ─────────────────────────────────────────────
# 'spectrum b' maps the B-factor values to a colour gradient.
# We use the official AlphaFold colour scheme:
#   orange (#FF7D45) → yellow (#FFDB13) → light blue (#65CBF3) → blue (#006EFF)
# This goes from low pLDDT (disordered) on the left to high (structured) right.
spectrum b, 0xFF7D45 0xFFDB13 0x65CBF3 0x006EFF, LOXL2, minimum=0, maximum=100

# ── Select and label known domains ───────────────────────────────────────────
# Signal peptide — disordered, low pLDDT
select signal_pep,  LOXL2 and resi 1-25

# SRCR domains 1 and 2 — structured
select SRCR_1,      LOXL2 and resi 91-175
select SRCR_2,      LOXL2 and resi 230-314

# Propeptide — disordered, keeps LOXL2 inactive before cleavage
select propeptide,  LOXL2 and resi 315-489

# LOX catalytic domain — very confident, the main drug target
select LOX_domain,  LOXL2 and resi 490-774

# Show propeptide as a dashed trace to indicate disorder
show  dash,  propeptide
hide  cartoon, propeptide

# Label key domain centroids
label LOXL2 and name CA and resi 130, "SRCR-1"
label LOXL2 and name CA and resi 272, "SRCR-2"
label LOXL2 and name CA and resi 400, "Propeptide\n(disordered)"
label LOXL2 and name CA and resi 630, "LOX catalytic"

set label_size,  12
set label_color, black

# ── Surface overlay ───────────────────────────────────────────────────────────
show surface, LOXL2
set transparency, 0.45, LOXL2

# ── Render ────────────────────────────────────────────────────────────────────
orient LOXL2
zoom   LOXL2, 3
ray 1400, 1400
png LOXL2_pLDDT_coloured.png, dpi=300

print "LOXL2 pLDDT visualisation complete! → LOXL2_pLDDT_coloured.png"
