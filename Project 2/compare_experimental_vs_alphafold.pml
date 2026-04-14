# compare_experimental_vs_alphafold.pml
# ─────────────────────────────────────────────────────────────────────────────
# Side-by-side comparison of the LOXL2 experimental crystal structure (5ZE3)
# and the AlphaFold prediction (Q9UBU1), after superimposition.
#
# Colour coding:
#   Marine blue   → experimental structure (the "ground truth")
#   Orange        → AlphaFold prediction
#   Red           → regions only in AlphaFold (propeptide + signal peptide)
#                   These are the regions crystallography cannot capture.
#
# Prerequisites:
#   - pdb5ze3.ent               (experimental — from 03_fetch_pdb_structures.py)
#   - LOXL2_AF_superimposed.pdb (rotated AF model — from 08_superimpose.py)
#
# How to run:  pymol -c compare_experimental_vs_alphafold.pml
# Output:      LOXL2_comparison.png
# ─────────────────────────────────────────────────────────────────────────────

# ── Load both structures ──────────────────────────────────────────────────────
# The experimental structure has only the LOX catalytic domain (partial).
load pdb5ze3.ent,               LOXL2_exp
load LOXL2_AF_superimposed.pdb, LOXL2_AF

# ── Display ───────────────────────────────────────────────────────────────────
hide everything
bg_color white
show cartoon, LOXL2_exp
show cartoon, LOXL2_AF
set cartoon_fancy_helices, 1

# ── Colours ───────────────────────────────────────────────────────────────────
# Experimental = marine blue (cold, solid, reliable)
# AlphaFold    = orange (warm, predicted)
color marine,  LOXL2_exp
color orange,  LOXL2_AF

# Highlight regions where AlphaFold has coverage but the crystal structure does not
# (these are the biologically important disordered regions that could not be crystallised)
select AF_extra, LOXL2_AF and resi 1-489
color red, AF_extra
# Show AF-only region as thicker trace to stand out
set cartoon_tube_radius, 0.4, AF_extra

# ── Positioning ───────────────────────────────────────────────────────────────
# Shift the two structures apart so both are fully visible
translate [-18, 0, 0], LOXL2_exp
translate [ 18, 0, 0], LOXL2_AF

# ── Labels ────────────────────────────────────────────────────────────────────
pseudoatom lbl_exp, pos=[-18, 35, 0]
pseudoatom lbl_af,  pos=[ 18, 35, 0]
pseudoatom lbl_red, pos=[ 18,  0, 0]

label lbl_exp, "Experimental\n(5ZE3 — LOX domain)"
label lbl_af,  "AlphaFold\n(Q9UBU1 — full length)"
label lbl_red, "Red = propeptide region\n(AlphaFold only)"

set label_size,  11
set label_color, black

# ── Render ────────────────────────────────────────────────────────────────────
# Wide panoramic format works well for side-by-side comparison
orient all
ray 2400, 900
png LOXL2_comparison.png, dpi=300

print "Comparison render complete! → LOXL2_comparison.png"
