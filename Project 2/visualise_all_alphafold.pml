# visualise_all_alphafold.pml
# ─────────────────────────────────────────────────────────────────────────────
# Loads all five AlphaFold fibrosis proteins side-by-side, coloured
# by pLDDT confidence score.
#
# This gives an at-a-glance comparison of which proteins are well-predicted
# (blue = structured) vs. which have large disordered regions (orange = flexible).
#
# Prerequisites:
#   Run 05_fetch_alphafold.py to download all five .pdb files.
#   All five files must be in the same directory as this script.
#
# How to run:  pymol -c visualise_all_alphafold.pml
# Output:      all_five_AF_plddt.png
# ─────────────────────────────────────────────────────────────────────────────

# ── Load all five AlphaFold models ───────────────────────────────────────────
load LOXL2_AF.pdb,  LOXL2
load CTHRC1_AF.pdb, CTHRC1
load ITGA11_AF.pdb, ITGA11
load SFRP2_AF.pdb,  SFRP2
load PDGFRB_AF.pdb, PDGFRB

# ── Display ───────────────────────────────────────────────────────────────────
hide everything
bg_color white
show cartoon, all
set cartoon_fancy_helices, 1
set cartoon_smooth_loops,  1

# ── Colour ALL by pLDDT (B-factor) ───────────────────────────────────────────
# Same colour scheme for all five — enables direct visual comparison
# of confidence levels across proteins
spectrum b, 0xFF7D45 0xFFDB13 0x65CBF3 0x006EFF, all, minimum=0, maximum=100

# ── Arrange in a row ─────────────────────────────────────────────────────────
# LOXL2 (774 aa) is the largest — give it more space
translate [  0, 0, 0], LOXL2
translate [ 60, 0, 0], CTHRC1
translate [110, 0, 0], ITGA11
translate [170, 0, 0], SFRP2
translate [230, 0, 0], PDGFRB

# ── Labels ────────────────────────────────────────────────────────────────────
label LOXL2  and name CA and resi 400, "LOXL2\n(774 aa)"
label CTHRC1 and name CA and resi 120, "CTHRC1\n(243 aa)"
label ITGA11 and name CA and resi 600, "ITGA11\n(1188 aa)"
label SFRP2  and name CA and resi 145, "SFRP2\n(295 aa)"
label PDGFRB and name CA and resi 550, "PDGFRB\n(1106 aa)"

set label_size,  13
set label_color, black

# ── Render ────────────────────────────────────────────────────────────────────
orient all
ray 3000, 700
png all_five_AF_plddt.png, dpi=300

print "All five AlphaFold models rendered! → all_five_AF_plddt.png"
