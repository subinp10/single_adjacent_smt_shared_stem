"""
routeJointNoAnimate.py

Standalone K/B/L extraction for a joint-model perTree (built via
jointMultiTargetLAFCA.buildPerTreeFromJointModel).

This is a direct lift of multiTargetLAFCA.multiTargetLAFCA()'s own
`animate=False` branch -- same parallel placement (placeOneTreeTimestep,
called once per active tree per timestep), same single combined DFL call
per timestep across ALL active trees, same post-routing relabel/wash
(relabelAfterRouting) -- just factored out so it can be called directly
on a perTree list without going through multiTargetLAFCA()'s own
tree-generation step (which the joint SMT solve replaces) and WITHOUT
the animate=True branch's matplotlib/GIF rendering.

No DP scheduling is involved anywhere in this path -- the timing that
drives the loop is whatever is already in each perTree entry's
timeStampMap (for the joint-model case, that's the joint solver's own
c_{k,i} cycle numbering, substituted in by buildPerTreeFromJointModel).
"""

import multiTargetLAFCA as mtl
from LAFCADFL.DFL import DFL as _realDFL


def routeJointNoAnimate(perTree):
    """
    perTree: list of dicts in multiTargetLAFCA's own shape (ratio, depth,
             M, assignment, parent, timeStampMap, ...) -- e.g. straight
             from jointMultiTargetLAFCA.buildPerTreeFromJointModel(...).

    Returns (K, B, L):
        K -- total flow count (number of DFL loading orders issued,
             i.e. total distinct reagent-loading events across all
             cycles/trees)
        B -- total bends across all loading paths
        L -- total path length across all loading paths
    """
    grid = [['*'] * mtl.GRID_COL for _ in range(mtl.GRID_ROW)]
    inlet, outlet = [0, mtl.GRID_ROW - 1], [mtl.GRID_COL - 1, 0]
    maxT = max(d["depth"] for d in perTree)
    totFlow = totBend = totLen = 0

    for t in range(maxT):
        combinedCellsToLoad = []
        combinedReagents = set()
        treeMixturesThisT = []  # (Mixtures, parentMix) pairs, relabeled AFTER routing

        for d in perTree:
            if t not in d["timeStampMap"]:
                continue
            mixList = d["timeStampMap"][t]
            Mixtures, loadingCells, reagentList, blockageList, units = \
                mtl.buildTimestepData(mixList, d["M"], d["assignment"])

            # placement only (LAFCA + corrected dfs blockage selection) --
            # no DFL call yet
            cellsToLoad, reagents = mtl.placeOneTreeTimestep(
                loadingCells, reagentList, blockageList, units, grid)
            combinedCellsToLoad.extend(cellsToLoad)
            combinedReagents.update(reagents)
            treeMixturesThisT.append((Mixtures, d["parent"]))

        # ---- ONE DFL call for this timestep, across ALL active trees ----
        if combinedCellsToLoad:
            loadingPaths = _realDFL(inlet, outlet, grid,
                                     list(combinedReagents), combinedCellsToLoad)
            for order in loadingPaths:
                totFlow += 1
                totBend += order[1]
                totLen += len(order[2])

        # ---- relabel/wash every active tree's box AFTER routing ----
        for Mixtures, parentMix in treeMixturesThisT:
            mtl.relabelAfterRouting(Mixtures, parentMix, grid)

    return totFlow, totBend, totLen
