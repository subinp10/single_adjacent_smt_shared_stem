"""
runSharedStem.py

Run this file directly to:
1. Generate the shared stem tree for two targets using the COMBINED
   single-SMT approach (stem + T1 + T2 solved together, R_c shared
   by construction).
2. Draw the tree PNG.
"""

from singleSMT import generateSharedStemTreeSingleSMT
from createTreeSharedStem import drawSharedStemTree

# -- configure here -----------------------------------------------------------
TARGET1   = [0.256, 0.316, 0.428]
TARGET2   = [0.316, 0.256, 0.428]
ERR       = 0.01
R         = 13
OUT_DIR   = "./sharedStem/"
LABEL     = "test"
MAX_DEPTH = 7
OUT_TREE  = "./sharedStem/sharedTree"
# -------------------------------------------------------------------------

print("Step 1: Generating shared stem tree (combined SMT)...")
result = generateSharedStemTreeSingleSMT(
    target1  = TARGET1,
    target2  = TARGET2,
    err      = ERR,
    outDir   = OUT_DIR,
    label    = LABEL,
    maxDepth = MAX_DEPTH,
)

if result is None:
    print("No shared stem found.")
elif "note" in result:
    print(result["note"])
else:
    c  = result['c']
    d1 = result['d1']
    d2 = result['d2']

    print(f"\nFound: c={c}  d1={d1}  d2={d2}")

    print(f"\nStep 2: Drawing tree -> {OUT_TREE}.png ...")
    drawSharedStemTree(result, OUT_TREE, R=R)
    print("Done.")
