# compare_all_proteins.pml
# ─────────────────────────────────────────────────────────────────────────────
# Loads all four cirrhosis-linked proteins side-by-side for visual comparison.
#
# Proteins displayed:
#   TGF-beta1 (1KTZ) — salmon
#   ACTA2     (3LUE) — blue
#   COL1A1    (1BKV) — orange (note: triple helix looks very different!)
#   TIMP1     (1UEA) — violet
#
# How to run:  pymol -c compare_all_proteins.pml
# Output:      four_proteins_comparison.png
# ─────────────────────────────────────────────────────────────────────────────

# ── Load all four structures ──────────────────────────────────────────────────
fetch 1KTZ, TGFb1
fetch 3LUE, ACTA2
fetch 1BKV, COL1A1
fetch 1UEA, TIMP1

# ── Clean display ─────────────────────────────────────────────────────────────
hide everything
show cartoon
bg_color white
set cartoon_fancy_helices, 1

# ── Unique colour per protein ─────────────────────────────────────────────────
# Different colours make the four proteins instantly distinguishable
color salmon,   TGFb1
color marine,   ACTA2
color orange,   COL1A1
color violet,   TIMP1

# ── Arrange side-by-side ──────────────────────────────────────────────────────
# 'translate' shifts each object along the X-axis by 30 Å increments
# so all four appear in a row without overlapping
translate [0,  0, 0], TGFb1
translate [30, 0, 0], ACTA2
translate [60, 0, 0], COL1A1
translate [90, 0, 0], TIMP1

# ── Labels ────────────────────────────────────────────────────────────────────
# Label the first CA atom of each object
# (these appear as floating text labels in the render)
label TGFb1  and name CA and resi 1, "TGF-beta1\n(1KTZ)"
label ACTA2  and name CA and resi 1, "ACTA2\n(3LUE)"
label COL1A1 and name CA and resi 1, "COL1A1\n(1BKV)"
label TIMP1  and name CA and resi 1, "TIMP1\n(1UEA)"

set label_size,  14
set label_color, black

# ── Render ────────────────────────────────────────────────────────────────────
# Wide panoramic format: 2400 wide × 600 tall
zoom all
ray 2400, 600
png four_proteins_comparison.png, dpi=300

# Save the PyMOL session (.pse) so you can reopen and rotate interactively
save cirrhosis_proteins_session.pse

print "Four-protein comparison complete! Saved → four_proteins_comparison.png"
