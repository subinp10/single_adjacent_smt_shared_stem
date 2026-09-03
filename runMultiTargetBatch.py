"""
runMultiTargetBatch.py

Headless batch runner for the independent multi-target pipeline (each
target gets its OWN Z3-solved tree via stg1's escalation, no cross-target
sharing objective) -- same underlying pipeline multiTargetGUI.py drives,
just without Tkinter and looped over every row of your input files
automatically.

IMPORTANT: this does NOT call multiTargetLAFCA.multiTargetLAFCA() as a
black box and then re-parse its output files afterward -- that was the
source of the earlier depth/fresh-reagent/union-cost bugs (guessed file
parsing that didn't match createTreeForShared's own conventions). Instead
it calls the SAME building blocks multiTargetLAFCA() itself calls
internally (generateSingleTree, cts.getMixerData, cts.getPlacement,
buildParentAndTimestamps), so:

    - depth is the ACTUAL depth generateSingleTree solved at (not
      re-derived from a text file afterward).
    - fresh reagent count is read directly from the x[i][r] array
      getMixerData returns, cross-referenced against each node's REAL
      solved cycle from timeStampMap -- "number of fresh reagent used at
      each timestamp, for each tree", summed across all trees/timestamps.
    - union_cost_independent is the number of DISTINCT (timestamp, reagent)
      pairs across ALL trees -- "number of different reagents used at
      each timestamp" -- again read from the real timeStampMap, not a
      guessed cycle formula.

K/B/L come from routeJointNoAnimate() -- the exact same parallel-DFL
routing loop multiTargetLAFCA()'s own animate=False branch runs -- applied
directly to the perTree list built here.

For each row i across the input files:
    1. Solves each target's own tree independently (own minimal depth,
       own z3 output file under OUT_ROOT/row<i+1>/tr<k>/).
    2. Renders each target's tree PNG (same as multiTargetLAFCA() does).
    3. Computes K, B, L via routeJointNoAnimate on the resulting perTree.
    4. Computes total_fresh_reagent and union_cost_independent directly
       from each tree's real x array + timeStampMap.
    5. Appends one row to the summary CSV, flushing after every row so
       the run can be resumed if interrupted.

Usage:
    python runMultiTargetBatch.py

Must be run from a directory where multiTargetLAFCA.py, stg1.py,
createTreeForShared.py, getLoadingData.py, routeJointNoAnimate.py, and
the LAFCADFL package are all importable.
"""

import os
import csv

import multiTargetLAFCA as mtl
import createTreeForShared as cts
from routeJointNoAnimate import routeJointNoAnimate

# =============================================================================
# CONFIG -- edit these before running
# =============================================================================

INPUT_FILES = [
    "t1_5.csv",
    "t2_5.csv",
    "t3_5.csv"
    

    
    
]  # row i across all of them = one batch of targets

ERR = 0.01    # error tolerance, applied to every target/row

OUT_ROOT = "./out2_5_fulldata"
SUMMARY_CSV = os.path.join(OUT_ROOT, "multiTargetBatchSummary.csv")


# =============================================================================
# Helpers
# =============================================================================

def readAllRows(path):
    with open(path, newline="") as f:
        r = csv.reader(f)
        next(r)  # header
        rows = []
        for row in r:
            if not row or all(v.strip() == "" for v in row):
                continue
            rows.append([float(x) for x in row])
        return rows


def loadCompletedRows(summaryPath):
    done = set()
    if os.path.exists(summaryPath):
        with open(summaryPath, newline="") as f:
            for row in csv.DictReader(f):
                done.add(int(row["row"]))
    return done


def buildIndependentPerTree(targets, rowDir):
    """
    Mirrors multiTargetLAFCA()'s own per-target loop (generateSingleTree ->
    getMixerData -> getPlacement -> buildParentAndTimestamps -> saveTree)
    exactly, but returns perTree (+ each tree's real x array) directly to
    the caller instead of routing/returning only aggregated K/B/L.

    Returns (perTree, xByTree, skipped):
        perTree  -- same shape multiTargetLAFCA()/routeJointNoAnimate expect
        xByTree  -- list of x[i][r] arrays, same order/index as perTree
        skipped  -- list of {"index", "ratio"} for targets with no SAT tree
    """
    perTree = []
    xByTree = []
    skipped = []

    for k, ratio in enumerate(targets):
        treeDir = os.path.join(rowDir, f"tr{k + 1}")
        z3opFile, depth = mtl.generateSingleTree(ratio, ERR, k, treeDir)
        if z3opFile is None:
            skipped.append({"index": k, "ratio": ratio})
            continue

        R = len(ratio)
        Rvals, M, x, W = cts.getMixerData(z3opFile, R, depth)
        assignment = cts.getPlacement(W, depth, list(mtl.SLOTS[k]))
        parent, timeStampMap = mtl.buildParentAndTimestamps(M, assignment, depth)

        treeImagePath = os.path.join(treeDir, f"tree{k + 1}")
        try:
            dot = cts.saveTree(x, W, treeImagePath, R=Rvals, N=mtl.N_MIXER)
            dot.render(treeImagePath, format='png', cleanup=True)
            treeImagePath = treeImagePath + ".png"
        except Exception as e:
            print(f"  [warn] tree diagram render failed for target {k+1}: {e}")
            treeImagePath = None

        perTree.append({
            "ratio": ratio, "depth": depth, "M": M,
            "assignment": assignment, "parent": parent,
            "timeStampMap": timeStampMap,
            "estK": 0, "estB": 0, "estL": 0,
            "treeImagePath": treeImagePath,
        })
        xByTree.append(x)

    return perTree, xByTree, skipped


def computeFreshAndUnionCost(perTree, xByTree):
    """
    Walks each tree's REAL solved timeStampMap (timestamp -> list of
    'M{i}' mix ids active at that timestamp) and cross-references each
    node's fresh-reagent load x[i][r] -- exactly "number of fresh reagent
    used at each timestamp, for each tree" (totalFresh), and "number of
    different reagents used at each timestamp" across ALL trees combined
    (union_cost_independent, counting each (timestamp, reagent) pair once
    regardless of how many trees/nodes load it at that same timestamp).
    """
    totalFresh = 0
    cycleReagentUsed = set()

    for tree, x in zip(perTree, xByTree):
        for t, mixList in tree["timeStampMap"].items():
            for mixId in mixList:
                if not mixId.startswith("M"):
                    continue
                try:
                    i = int(mixId[1:])
                except ValueError:
                    continue
                if i >= len(x):
                    continue  # not a fresh-reagent-bearing node index
                for r, amt in enumerate(x[i]):
                    if amt > 0:
                        totalFresh += amt
                        cycleReagentUsed.add((t, r))

    return totalFresh, len(cycleReagentUsed)


# =============================================================================
# Main
# =============================================================================

def main():
    os.makedirs(OUT_ROOT, exist_ok=True)

    perFileRows = [readAllRows(p) for p in INPUT_FILES]
    print("Input file row counts:")
    for p, rows in zip(INPUT_FILES, perFileRows):
        print(f"  {p}: {len(rows)} rows")

    nRows = min(len(rows) for rows in perFileRows)
    if any(len(rows) != nRows for rows in perFileRows):
        print(f"[warn] input files have differing row counts; "
              f"truncating to shortest = {nRows}")

    K = len(INPUT_FILES)
    fieldnames = (["row", "status", "K", "B", "L", "total_fresh_reagent",
                   "union_cost_independent"] +
                  [f"depth_t{k+1}" for k in range(K)] + ["skipped"])

    alreadyDone = loadCompletedRows(SUMMARY_CSV)
    writeHeader = not os.path.exists(SUMMARY_CSV)

    with open(SUMMARY_CSV, "a", newline="") as sf:
        writer = csv.DictWriter(sf, fieldnames=fieldnames)
        if writeHeader:
            writer.writeheader()
            sf.flush()

        for i in range(nRows):
            if i in alreadyDone:
                continue

            targets = [perFileRows[f][i] for f in range(K)]
            print(f"\n=== row {i + 1}/{nRows} === targets={targets}")

            rowDir = os.path.join(OUT_ROOT, f"row{i + 1}")
            os.makedirs(rowDir, exist_ok=True)

            perTree, xByTree, skipped = buildIndependentPerTree(targets, rowDir)

            record = {"row": i + 1}
            for k in range(K):
                record[f"depth_t{k+1}"] = ""

            if not perTree:
                record["status"] = "no_sat_tree_any_target"
                record["K"] = record["B"] = record["L"] = ""
                record["total_fresh_reagent"] = ""
                record["union_cost_independent"] = ""
                record["skipped"] = "all"
                writer.writerow(record)
                sf.flush()
                continue

            record["status"] = "ok"
            record["skipped"] = "; ".join(str(s["index"] + 1) for s in skipped) or ""

            # depth per target, in ORIGINAL target order (skipped ones stay blank)
            solvedIdx = [k for k in range(K)
                         if k not in {s["index"] for s in skipped}]
            for pos, k in enumerate(solvedIdx):
                record[f"depth_t{k+1}"] = perTree[pos]["depth"]

            try:
                routeK, routeB, routeL = routeJointNoAnimate(perTree)
                record["K"], record["B"], record["L"] = routeK, routeB, routeL
            except Exception as e:
                print(f"  [warn] routing failed for row {i+1}: {e}")
                record["K"] = record["B"] = record["L"] = ""

            totalFresh, unionCost = computeFreshAndUnionCost(perTree, xByTree)
            record["total_fresh_reagent"] = totalFresh
            record["union_cost_independent"] = unionCost

            writer.writerow(record)
            sf.flush()

    print(f"\nDone. Summary: {SUMMARY_CSV}")


if __name__ == "__main__":
    main()