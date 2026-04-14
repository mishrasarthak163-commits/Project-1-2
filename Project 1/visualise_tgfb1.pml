# visualise_tgfb1.pml
# ─────────────────────────────────────────────────────────────────────────────
# Visualises TGF-beta1 (PDB: 1KTZ) coloured by secondary structure.
#
# How to run:
#   Option A — from PyMOL GUI console:  run visualise_tgfb1.pml
#   Option B — from terminal:           pymol -c visualise_tgfb1.pml
#   Option C — headless (no GUI):       pymol -cq visualise_tgfb1.pml
#
# Output: TGF_beta1_structure.png (1200×1200, ray-traced)
# ─────────────────────────────────────────────────────────────────────────────

# ── Step 1: Fetch the structure directly from RCSB ───────────────────────────
# 'fetch' downloads the PDB file and loads it in one command.
# The second argument ('TGF_beta1') is the name we give the object in PyMOL.
fetch 1KTZ, TGF_beta1

# ── Step 2: Set a clean starting display ─────────────────────────────────────
hide everything, TGF_beta1      # start from blank — no default lines
bg_color white                  # white background is clearest for publications
set ray_opaque_background, on   # solid white (not transparent) when saving

# ── Step 3: Show as cartoon ───────────────────────────────────────────────────
# 'cartoon' renders helices as ribbons/tubes and sheets as arrows.
# It is the most informative representation for seeing overall fold.
show cartoon, TGF_beta1
set cartoon_fancy_helices, 1    # nicer helices with a lip/ridge
set cartoon_smooth_loops, 1     # smoother loop regions

# ── Step 4: Colour by secondary structure ────────────────────────────────────
#   ss h  = alpha helices  → red (warm colour = structured)
#   ss s  = beta strands   → yellow (sheets)
#   ss l  = loops / coils  → forest green (flexible regions)
color red,          TGF_beta1 and ss h
color yellow,       TGF_beta1 and ss s
color forest,       TGF_beta1 and ss l+''

# ── Step 5: Highlight the receptor-binding loop ───────────────────────────────
# This loop (residues 50–100) is where TGF-beta1 contacts its receptor.
# It is the primary drug target in fibrosis therapy.
select binding_loop, TGF_beta1 and resi 50-100
color cyan, binding_loop
show sticks, binding_loop       # show individual bonds — useful for drug design

# ── Step 6: Add a semi-transparent surface ───────────────────────────────────
# The surface shows the molecular shape — what a drug or receptor would "see".
show surface, TGF_beta1
set transparency, 0.5, TGF_beta1    # 50% transparent so cartoon is still visible

# ── Step 7: Final framing and render ─────────────────────────────────────────
orient TGF_beta1       # auto-orient to best viewing angle
zoom TGF_beta1, 5      # zoom out slightly so the whole protein is visible

# ray-trace produces publication-quality anti-aliased images
# 1200x1200 pixels is standard for journal figures
ray 1200, 1200
png TGF_beta1_structure.png, dpi=300

print "TGF-beta1 visualisation complete! Saved → TGF_beta1_structure.png"
