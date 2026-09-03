"""
multiTargetSkewed.py

Multi-target skewed sample preparation via *repeated single-tree generation*.

Rather than solving all target ratios in one monolithic SMT model, this driver
generates an independent skewed mixing tree for each target by repeatedly
invoking the existing single-tree Z3 pipeline (skewedTreeGenerator), then
assembles the trees onto one chip and reports combined metrics.

Per run:
    targets --> [ generateSingleTree() x K ]      (repeated single-tree gen)
            --> reconstruct each tree
            --> static placement on a 10x10 grid
            --> aggregate area / waste / reagent usage / loading cost (K, B, L)
            --> optional graphical simulation (GIF) of the parallel loading

Reuses the project's existing modules (must be on the import path):
    skewedTreeGenerator : SMT model emission (initOPT, addvariables, ...)
    createTreeForShared : getMixerData, getPlacement
    getLoadingData      : getLoadingData (K, B, L)   [pulls NTM via createTree]
    LAFCADFL            : DFL  (parallel loading + the graphical simulation)

Requirements: z3-solver, and the NTM / LAFCADFL modules used elsewhere in the
project. graphviz is NOT required unless render=True (tree.png per target).
matplotlib is only required if animate=True.

Windows-safe: launches the generated solver with sys.executable, never "python3".

NOTE on placement compatibility: getPlacement's column-shift heuristic for
non-adjacent sharing is kept EXACTLY as originally written (not the
interval-based holding-cell redesign explored elsewhere), because
getLoadingData's K/B/L lookup tables (KBL_0_*, KBL_1_*) are keyed to the
EXACT three cell offsets that shift scheme produces. Changing getPlacement
would silently make getLoadingData return wrong numbers without erroring.
"""

import os
import sys
import csv
import subprocess
from LAFCADFL.DFL import DFL
import stg1 as stg
import createTreeForShared as cts
import getLoadingData as gld


# ----- chip / mixer configuration -------------------------------------------
GRID_ROW, GRID_COL = 10, 10
N_MIXER = 4  # mixer-4

# Up to 9 static 2x2 mixer slots on a 10x10 grid (same convention as multiTarget.py)
SLOTS = [[1, 1], [4, 1], [4, 4], [1, 4], [1, 7], [4, 7], [7, 1], [7, 4], [7, 7]]
# 2x2 block offsets for the four cells of one mixer
BLOCK_X = [0, 0, 1, 1]
BLOCK_Y = [0, 1, 1, 0]


# ----- 1. repeated single-tree generation -----------------------------------
def generateSingleTree(targetRatio, err, ind, directory, N=N_MIXER, maxDepth=10):
    """
    Emit + solve ONE single-target skewed-tree SMT instance, escalating the tree
    depth from 2 upward until SAT (or until maxDepth is reached).

    This is the "repeated single tree generation" core: call it once per target.

    Returns
        (z3opFile, depth)  on success
        (None, None)       if no SAT tree exists at depth < maxDepth
    """
    os.makedirs(directory, exist_ok=True)

    R = len(targetRatio)
    # forward slashes: this path is embedded verbatim into generated Python source
    z3opFile = (directory + f"/z3outputFile_{ind}_{err}.txt").replace("\\", "/")
    scriptFile = (directory + f"/_genZ3_{ind}.py").replace("\\", "/")

    d = 2
    while d < maxDepth:
        # generate the SMT program for the current depth
        opfile = open(scriptFile, "w+")
        stg.initOPT(opfile)
        stg.addvariables(d, R, opfile)
        stg.linearityVariables(d, N, R, opfile)
        stg.nonNegativityConstraints(d, N, R, opfile)
        stg.mixerConsistencyConstraints(d, N, R, opfile)
        stg.fluidSharingConstraint(d, N, R, opfile)
        stg.setTarget(targetRatio, err, R, N, d, opfile)
        stg.finishOPT(opfile, z3opFile)   # finishOPT closes opfile

        # solve it with the SAME interpreter (Windows-safe)
        subprocess.call([sys.executable, scriptFile])

        with open(z3opFile, "r") as f:
            content = f.read().strip()

        if content == "unsat":
            d += 1
        else:
            return z3opFile, d

    return None, None


# ----- 2. per-tree reconstruction & metrics ---------------------------------
def boundingBoxArea(cells):
    """Rectangular area spanned by a list of [row, col] cells."""
    if not cells:
        return 0
    xs = [c[0] for c in cells]
    ys = [c[1] for c in cells]
    return (max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)


def computeProtectedCells(W, assignment, depth):
    """
    For every non-adjacent share (i > j+1), identify which of mixer i's
    OWN getPlacement-assigned cells hold the droplet destined for mixer j,
    and the (start, end) timestep window they must stay protected
    (untouched by anyone else) between i producing it and j consuming it.

    Uses the EXISTING placement algorithm as-is (no new cell geometry, no
    separate holding area) -- just tracks, per cell, which timestep window
    nothing is allowed to overwrite it in. Blocking only applies STRICTLY
    BETWEEN start and end: the source mixer's own write at t == start (its
    own production) and the consumer's write at t == end (its own genuine
    consumption) both always go through -- only an unrelated intervening
    mixer's write, at start < t < end, gets blocked.

    Returns {(r, c): (startTimestep, endTimestep)}.
    """
    protected = {}
    for i in W:
        for j in W[i]:
            if i > j + 1:
                units = W[i][j]
                start = depth - 1 - i
                end = depth - 1 - j
                sourceCells = assignment[f"M{i}"][1]
                # same priority LAFCA_DFL uses for partial holds: bottom-left
                # first, then top-left, then top-right
                priority = [sourceCells[2], sourceCells[0], sourceCells[1]]
                for cell in priority[:units]:
                    key = tuple(cell)
                    # if a cell is claimed by more than one pending hold,
                    # protect it through whichever consumer runs latest
                    cur = protected.get(key)
                    if cur is None or end > cur[1]:
                        protected[key] = (start, end)
    return protected


def treeMetrics(z3opFile, depth, R, slot):
    """
    Reconstruct a solved tree and compute its metrics, placed at `slot`.
    K/B/L here comes straight from getLoadingData -- the real, validated
    lookup-table cost estimator -- NOT from re-running DFL live.
    No graphviz/dot dependency (waste/reagent computed directly from x, W).
    """
    _, M, x, W = cts.getMixerData(z3opFile, R, depth)

    # waste = unused sharing capacity at nodes that share fluid
    # (same definition as createTreeForShared.saveTree_getArea)
    waste = 0
    for i in W:
        waste += N_MIXER - sum(W[i].values())

    # total units of pure reagent loaded across all mixers
    reagentUsage = sum(sum(node) for node in x)

    # placement of this tree at its slot -> occupied cells
    assignment = cts.getPlacement(W, depth, list(slot))
    cells = []
    for mix in assignment:
        cells.extend(assignment[mix][1])

    # which of getPlacement's own assigned cells must be protected from
    # being overwritten before their non-adjacent consumer's timestep
    protected = computeProtectedCells(W, assignment, depth)

    # loading cost (K, B, L) for this tree at its slot -- the real
    # getLoadingData lookup-table estimator
    K, B, L = gld.getLoadingData(
        z3opFile, list(slot), depth, GRID_ROW, GRID_COL, N_MIXER, R
    )

    # per-tree footprint: same two numbers generateSkewedTree()'s single-target
    # path gets from saveTree_getArea's ld.boundingbox() call -- BoundingBox
    # (count of unique occupied cells) and area (rectangular bounding-box
    # area). This is THIS tree's own footprint, not the batch-combined one
    # computed later in prepareMultiTarget from every tree's cells together.
    boundingBox = len({tuple(c) for c in cells})
    area = boundingBoxArea(cells)

    return {
        "waste": waste,
        "reagent": reagentUsage,
        "cells": cells,
        "boundingBox": boundingBox, "area": area,
        "K": K, "B": B, "L": L,
        "M": M, "x": x, "W": W, "assignment": assignment, "protected": protected,
    }


# ----- 3. parallel loading (also feeds the graphical simulation) ------------
def buildTimelines(perTree):
    """
    Build, per tree, a timeline indexed by GLOBAL timestep t, using each
    mixer's ACTUAL cells from cts.getPlacement (the existing placement
    algorithm, unchanged, shift included). Non-adjacent sharing is handled
    by BLOCKING overwrites to protected cells (see computeProtectedCells /
    treeMetrics), not by relocating anything -- this uses the existing
    placement as-is, no new cell geometry.

    Returns (timelines, maxT) where timelines is a list of
    (slotIdx, timeline), and timeline[t] is either None (no mixer active
    for this tree at t) or {"cells": [...], "droplets": [...]}.
    """
    timelines = []
    maxT = 0
    for (slotIdx, _ratio, depth, m) in perTree:
        M = m["M"]
        assignment = m["assignment"]  # {f"M{i}": [timestamp, cells]}
        timeline = [None] * depth
        for mixerName, (ts, cells) in assignment.items():
            timeline[ts] = {"cells": [c[:] for c in cells],
                             "droplets": M.get(mixerName, [])}
        timelines.append((slotIdx, timeline))
        maxT = max(maxT, depth)
    return timelines, maxT


def parallelLoading(perTree):
    """
    Route droplet loading for all trees *simultaneously*, timestamp by timestamp,
    using the project's DFL router. Returns (flow, bendings, pathLength, collisions)
    or None if the router is unavailable.

    Non-adjacent sharing: uses getPlacement's EXISTING cell assignment as-is
    (no new cell geometry). A cell holding a droplet destined for a future
    non-adjacent consumer is PROTECTED (write-blocked) until that consumer's
    own timestep actually arrives -- clear it away only once it's genuinely
    no longer needed by anyone; never wash a cell mid-use for an adjacent
    share, since that's already just the natural next-timestep overwrite.

    Cross-tree collision handling: since getPlacement's shift for non-adjacent
    sharing can push a mixer's cells outside its tree's default SLOTS footprint,
    two DIFFERENT trees can end up claiming the same cell in the same timestep.
    Rather than silently letting the second tree's write clobber the first
    (which was the previous behavior -- undetectable data loss), the second
    write is BLOCKED and the collision is counted and logged.
    """
    try:
        from LAFCADFL import DFL
    except Exception as e:
        print(f"  (parallel loading skipped: {e})")
        return None

    grid = [['*'] * GRID_COL for _ in range(GRID_ROW)]
    timelines, maxT = buildTimelines(perTree)
    protectedBySlot = {slotIdx: dict(m["protected"]) for (slotIdx, _r, _d, m) in perTree}

    totFlow = totBend = totLen = totCollisions = totBlocked = 0
    for t in range(maxT):
        cellsToLoad = []
        required = set()
        claimedBy = {}   # (r,c) -> slotIdx that already placed a droplet here THIS timestep
        for slotIdx, timeline in timelines:
            if t >= len(timeline) or timeline[t] is None:
                continue
            entry = timeline[t]
            cells, droplets = entry["cells"], entry["droplets"]
            protectedUntil = protectedBySlot[slotIdx]
            for j, elm in enumerate(droplets[:4]):
                r, c = cells[j]
                key = (r, c)

                # --- protect cells still needed by a FUTURE non-adjacent
                # consumer: block only if this write is STRICTLY between the
                # hold's start and end -- the source mixer's own production
                # (t==start) and the consumer's own consumption (t==end)
                # always go through ---
                rng = protectedUntil.get(key)
                if rng is not None:
                    start, endTs = rng
                    if start < t < endTs:
                        totBlocked += 1
                        print(f"  [protected] t={t} slot{slotIdx+1} cell={key}: still needed "
                              f"for non-adjacent sharing until t={endTs} -- '{elm}' write blocked")
                        continue
                    if t == endTs:
                        del protectedUntil[key]   # consumed now -- free for reuse from here on

                if key in claimedBy and claimedBy[key] != slotIdx:
                    totCollisions += 1
                    print(f"  [collision] t={t} cell={key}: slot{claimedBy[key]+1} already "
                          f"occupies this cell -- slot{slotIdx+1}'s '{elm}' was blocked")
                    continue  # do NOT overwrite the cell that got there first
                claimedBy[key] = slotIdx
                grid[r][c] = elm
                if isinstance(elm, str) and elm.startswith('R'):
                    cellsToLoad.append([r, c])
                    required.add(elm)

        if not cellsToLoad:
            continue
        reagents = list(required)
        paths = DFL.DFL([0, GRID_ROW - 1], [GRID_COL - 1, 0], grid, reagents, cellsToLoad)
        for order in paths:
            totFlow += 1
            totBend += order[1]
            totLen += len(order[2])

    if totCollisions:
        print(f"  ** {totCollisions} cross-tree cell collision(s) blocked this batch **")
    if totBlocked:
        print(f"  ** {totBlocked} write(s) blocked to protect pending non-adjacent shares **")
    return totFlow, totBend, totLen, totCollisions


# ----- 3b. graphical simulation (replaces the missing repeatskew_simulation) -
_REAGENT_PALETTE = ["#D85A30", "#378ADD", "#639922", "#D4537E",
                    "#EF9F27", "#7F77DD", "#1D9E75", "#F09595"]


def _reagentColor(label):
    try:
        idx = int(label[1:]) - 1
    except (ValueError, IndexError):
        idx = 0
    return _REAGENT_PALETTE[idx % len(_REAGENT_PALETTE)]


def animate_batch(perTree, outfile="batch_sim.gif", fps=2,
                   grid_row=GRID_ROW, grid_col=GRID_COL):
    """
    Graphical simulation of the whole batch's parallel loading, timestep by
    timestep, on the shared 10x10 chip. Built directly on the SAME grid/
    timeline construction parallelLoading() already uses, so the animation
    and the reported parallel K/B/L numbers are guaranteed consistent with
    each other.

    perTree : the list returned by prepareMultiTarget's per-target loop --
              [(slotIdx, ratio, depth, metricsDict), ...]
    """
    from LAFCADFL import DFL
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import matplotlib.animation as animation
    from matplotlib.animation import PillowWriter

    INK, MUTE = "#2b2f36", "#5b6472"
    GRID_DOT, EMPTY = "#d7dce4", "#f2f1ec"
    MIXER_COLOR, PATH_GLOW = "#F2A623", "#ffffff"
    TARGET_COLOR = "#7F4FDB"  # distinct color for the final root (M0) mixed result
    BOX_COLORS = ["#2b2f36", "#993c1d", "#3c3489", "#085041",
                  "#712b13", "#26215c", "#04342c", "#4a1b0c", "#4b1528"]
    INLET, OUTLET = (0, grid_row - 1), (grid_col - 1, 0)

    timelines, maxT = buildTimelines(perTree)
    ratioBySlot = {slotIdx: ratio for (slotIdx, ratio, _d, _m) in perTree}
    depthBySlot = {slotIdx: depth for (slotIdx, _r, depth, _m) in perTree}
    protectedBySlot = {slotIdx: dict(m["protected"]) for (slotIdx, _r, _d, m) in perTree}

    # ----- build frames: one incremental reveal per routed droplet -----
    frames = []
    grid = [['*'] * grid_col for _ in range(grid_row)]      # fed to DFL -- may
                                                               # be mutated by DFL
                                                               # as a side effect
                                                               # of routing
    visGrid = [['*'] * grid_col for _ in range(grid_row)]   # what actually gets
                                                               # rendered -- NEVER
                                                               # read from `grid`
                                                               # after the DFL call,
                                                               # only written to
                                                               # explicitly below
    cumFlow = cumBend = cumLen = 0
    totCollisions = 0

    for t in range(maxT):
        cellsToLoad, required = [], set()
        cellLabel = {}   # (r,c) -> reagent label expected there this timestep
        activeSlots = []
        activeBoxes = {}  # slotIdx -> (boxCells, stampLabel) for the mix-stamp step
        protectedCellsThisFrame = []  # currently-protected cells, for the outline marker
        claimedBy = {}    # (r,c) -> slotIdx that already placed a droplet here THIS timestep
        collisions = []   # cross-tree collisions detected this timestep, for the frame marker

        for slotIdx, timeline in timelines:
            if t >= len(timeline) or timeline[t] is None:
                continue
            activeSlots.append(slotIdx)
            depth = len(timeline)
            mixerIdx = depth - 1 - t   # t=0 is the deepest/leaf mixer
            entry = timeline[t]
            cells, droplets = entry["cells"], entry["droplets"]
            boxCells = [(r, c) for r, c in cells]
            protectedUntil = protectedBySlot[slotIdx]
            # "T" marks the FINAL result (mixerIdx==0 is the root, M0 --
            # the target ratio has been reached); otherwise this timestep's
            # mixed output becomes intermediate fluid "M{mixerIdx}"
            activeBoxes[slotIdx] = (boxCells, "T" if mixerIdx == 0 else f"M{mixerIdx}")

            for j, elm in enumerate(droplets[:4]):
                r, c = cells[j]
                key = (r, c)

                # --- protect cells still needed by a FUTURE non-adjacent
                # consumer: block only if strictly between the hold's start
                # and end -- the source mixer's own production (t==start)
                # and the consumer's own consumption (t==end) always go
                # through ---
                rng = protectedUntil.get(key)
                if rng is not None:
                    start, endTs = rng
                    if start < t < endTs:
                        print(f"  [protected] t={t} slot{slotIdx+1} cell={key}: still needed "
                              f"for non-adjacent sharing until t={endTs} -- '{elm}' write blocked")
                        continue
                    if t == endTs:
                        del protectedUntil[key]   # consumed now -- free for reuse from here on

                # --- cross-tree collision check: block, don't silently overwrite ---
                if key in claimedBy and claimedBy[key] != slotIdx:
                    ownerSlot = claimedBy[key]
                    totCollisions += 1
                    collisions.append({"cell": key, "ownerSlot": ownerSlot,
                                        "blockedSlot": slotIdx, "elm": elm})
                    print(f"  [collision] t={t} cell={key}: slot{ownerSlot+1} already "
                          f"occupies this cell -- slot{slotIdx+1}'s '{elm}' was blocked")
                    continue  # do NOT overwrite the cell that got there first

                claimedBy[key] = slotIdx
                grid[r][c] = elm
                if isinstance(elm, str) and elm.startswith('R'):
                    cellsToLoad.append([r, c])
                    required.add(elm)
                    cellLabel[(r, c)] = elm
                elif isinstance(elm, str) and elm.startswith('M'):
                    visGrid[r][c] = elm   # carried droplet: not routed, appears instantly

            # any cell STILL protected right now (start < t < end) gets
            # flagged for the outline marker
            for key, (start, endTs) in protectedUntil.items():
                if start < t < endTs and key in [tuple(c) for c in boxCells]:
                    protectedCellsThisFrame.append({"slotIdx": slotIdx, "cell": key})

        treeStatus = [{"idx": s, "active": s in activeSlots, "done": t >= len(tl),
                        "box": activeBoxes.get(s, (None, None))[0],
                        "protectedCells": [pc["cell"] for pc in protectedCellsThisFrame
                                            if pc["slotIdx"] == s]}
                       for s, tl in timelines]

        if cellsToLoad:
            reagents = list(required)
            paths = DFL.DFL([0, grid_row - 1], [grid_col - 1, 0], grid, reagents, cellsToLoad)

            for order in paths:
                label, bends, path = order[0], order[1], order[2]
                cumFlow += 1
                cumBend += bends
                cumLen += len(path)
                # show the full path in 3 sub-frames: 33%, 66%, 100%
                # This guarantees the outlet is always reached before
                # advancing to the next flow, regardless of path length
                step = max(1, len(path) // 3)
                for kk in range(step, len(path) + 1, step):
                    headPath = path[:kk]
                    for cell in headPath:
                        tc = (cell[0], cell[1])
                        if cellLabel.get(tc) == label:
                            visGrid[tc[0]][tc[1]] = label
                    frames.append({
                        "grid": [row[:] for row in visGrid],
                        "paths": [(label, headPath)],
                        "t": t, "treeStatus": treeStatus, "collisions": collisions,
                        "cum": (cumFlow, cumBend, cumLen),
                    })
                # ensure destination is marked
                for cell in path:
                    tc = (cell[0], cell[1])
                    if cellLabel.get(tc) == label:
                        visGrid[tc[0]][tc[1]] = label
                # hold one extra frame showing the COMPLETE path inlet→outlet
                # so the full route is always clearly visible before advancing
                frames.append({
                    "grid": [row[:] for row in visGrid],
                    "paths": [(label, path)],   # full path, always reaches outlet
                    "t": t, "treeStatus": treeStatus, "collisions": collisions,
                    "cum": (cumFlow, cumBend, cumLen),
                })

        # ----- mixing complete for this timestep: STAMP the result --------
        # Every active mixer's 4 cells collapse into ONE uniform label --
        # replacing whatever individual reagent letters / carried tokens
        # were there. This is the actual physical mixing action: 4 separate
        # droplets combine into a single homogeneous output. Do NOT leave
        # raw reagent letters showing after this point -- only the mixer's
        # own intermediate label (or "T" once the root, M0, is reached)
        # should ever be visible from here on for these cells.
        # EXCEPTION: never stamp over a cell still protected for a future
        # non-adjacent consumer -- that content isn't this mixer's to claim.
        for slotIdx, (boxCells, stampLabel) in activeBoxes.items():
            protectedUntil = protectedBySlot[slotIdx]
            for (r, c) in boxCells:
                rng = protectedUntil.get((r, c))
                if rng is not None and rng[0] < t < rng[1]:
                    continue
                visGrid[r][c] = stampLabel

        # hold a few frames on the mixed/stamped result before advancing
        for _ in range(3):
            frames.append({"grid": [row[:] for row in visGrid], "paths": [], "t": t,
                            "treeStatus": treeStatus, "collisions": collisions,
                            "cum": (cumFlow, cumBend, cumLen)})

        # ----- clear cells not needed in any future timestep ---------------
        # A cell can be cleared if:
        #   1. Its mixer's timeline is now done (no more timesteps for this tree)
        #   2. It is NOT protected for a future non-adjacent share
        # This matches the physical reality: once mixing is complete and the
        # output has been either consumed or held, all other cells in that
        # mixer's box are washed clear and available for future use.
        for slotIdx, tl in timelines:
            if t >= len(tl) - 1:   # this tree's last active timestep just finished
                base = SLOTS[slotIdx]
                boxCells = [(base[0] + BLOCK_X[j], base[1] + BLOCK_Y[j]) for j in range(4)]
                protectedUntil = protectedBySlot[slotIdx]
                for (r, c) in boxCells:
                    key = (r, c)
                    rng = protectedUntil.get(key)
                    if rng is not None:
                        continue   # still needed for non-adjacent share -- keep it
                    visGrid[r][c] = '*'   # wash clear
                    grid[r][c] = '*'
            elif t < len(tl) - 1:
                # tree still has future timesteps -- only clear cells whose
                # mixer is done AND that aren't protected AND won't be reused
                # by this tree's own next timestep (the box is reused each
                # timestep in a skewed tree, so we only clear if nothing is
                # held for a non-adjacent share)
                base = SLOTS[slotIdx]
                boxCells = [(base[0] + BLOCK_X[j], base[1] + BLOCK_Y[j]) for j in range(4)]
                protectedUntil = protectedBySlot[slotIdx]
                for (r, c) in boxCells:
                    key = (r, c)
                    rng = protectedUntil.get(key)
                    if rng is not None and t < rng[1]:
                        continue   # protected until a future timestep -- keep it
                    # not protected -- will be overwritten next timestep anyway
                    # but clear visGrid now so it looks clean between timesteps
                    visGrid[r][c] = '*'
                    grid[r][c] = '*'

    # ----- render -----
    def xy(r, c):
        return c, (grid_row - 1 - r)

    def cellColor(elm):
        if not isinstance(elm, str) or elm == "*":
            return EMPTY, ""
        if elm == "T":
            return TARGET_COLOR, "T"
        if elm.startswith("R"):
            return _reagentColor(elm), elm
        if elm.startswith("M"):
            return MIXER_COLOR, elm
        return EMPTY, ""

    fig = plt.figure(figsize=(11, 6.4))
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.1, 1.1], figure=fig, wspace=0.12)
    ag = fig.add_subplot(gs[0, 0])
    asx = fig.add_subplot(gs[0, 1])

    def drawGrid(fr):
        ag.clear()
        ag.set_xlim(-0.7, grid_col - 0.3)
        ag.set_ylim(-0.7, grid_row - 0.3)
        ag.set_aspect("equal")
        ag.axis("off")
        ag.text((grid_col - 1) / 2, grid_row - 0.2, f"batch chip -- t = {fr['t']}",
                 ha="center", fontsize=11, fontweight="bold", color=INK)

        for r in range(grid_row):
            for c in range(grid_col):
                px, py = xy(r, c)
                ag.add_patch(mpatches.Circle((px, py), 0.07, facecolor=GRID_DOT,
                                              edgecolor="none", zorder=1))

        for r in range(grid_row):
            for c in range(grid_col):
                color, label = cellColor(fr["grid"][r][c])
                if not label:
                    continue
                px, py = xy(r, c)
                ag.add_patch(mpatches.FancyBboxPatch(
                    (px - 0.42, py - 0.42), 0.84, 0.84,
                    boxstyle="round,pad=0.02,rounding_size=0.12",
                    facecolor=color, edgecolor="white", linewidth=1.2, zorder=2))
                ag.text(px, py, label, ha="center", va="center", fontsize=7,
                        fontweight="bold", color="white", zorder=3)

        for label, path in fr["paths"]:
            if len(path) < 2:
                continue
            color = _reagentColor(label)
            pts = [xy(p[0], p[1]) for p in path]
            ag.plot([p[0] for p in pts], [p[1] for p in pts], color=color,
                     lw=3.2, alpha=0.85, zorder=4, solid_capstyle="round")
            hx, hy = pts[-1]
            ag.add_patch(mpatches.Circle((hx, hy), 0.2, facecolor=PATH_GLOW,
                                          edgecolor=color, lw=1.6, zorder=5))

        for st in fr["treeStatus"]:
            i = st["idx"]
            boxColor = BOX_COLORS[i % len(BOX_COLORS)]
            if st["box"] is not None:
                box = st["box"]   # real (possibly shifted) cells for the active mixer
            else:
                # idle/done: no active mixer this frame -- show a faint outline
                # at the tree's default home slot, purely as a visual anchor
                base = SLOTS[i]
                box = [(base[0] + BLOCK_X[j], base[1] + BLOCK_Y[j]) for j in range(4)]
            rs = [c[0] for c in box]
            cs = [c[1] for c in box]
            x0, y0 = xy(max(rs), min(cs))
            ag.add_patch(plt.Rectangle((x0 - 0.46, y0 - 0.46),
                         (max(cs) - min(cs)) + 0.92, (max(rs) - min(rs)) + 0.92,
                         fill=False, edgecolor=boxColor, lw=2.0, zorder=6,
                         alpha=1.0 if st["active"] else 0.25))
            ag.text(x0 - 0.46, y0 + 0.55, f"slot{i+1}", fontsize=7,
                    fontweight="bold", color=boxColor, zorder=6)

            for (pr, pc) in st.get("protectedCells", []):
                ppx, ppy = xy(pr, pc)
                ag.add_patch(plt.Rectangle((ppx - 0.46, ppy - 0.46), 0.92, 0.92,
                             fill=False, edgecolor=boxColor, lw=1.6, ls=(0, (4, 2)), zorder=7))

        # cross-tree collisions this timestep: a droplet was blocked from
        # writing to a cell another tree already claimed -- shown as a red
        # warning ring so it's never silently invisible
        for coll in fr.get("collisions", []):
            r, c = coll["cell"]
            px, py = xy(r, c)
            ag.add_patch(mpatches.Circle((px, py), 0.46, fill=False,
                                          edgecolor="#E24B4A", lw=2.4, zorder=8))
            ag.text(px, py - 0.62, f"blocked: slot{coll['blockedSlot']+1}",
                    ha="center", fontsize=6.5, fontweight="bold",
                    color="#E24B4A", zorder=8)

        ix, iy = xy(*INLET)
        ox, oy = xy(*OUTLET)
        ag.add_patch(mpatches.Circle((ix, iy), 0.28, facecolor=INK, edgecolor="none", zorder=5))
        ag.text(ix, iy + 0.55, "INLET", ha="center", fontsize=7.5, fontweight="bold", color=INK)
        ag.add_patch(mpatches.Circle((ox, oy), 0.28, facecolor=INK, edgecolor="none", zorder=5))
        ag.text(ox, oy - 0.55, "OUTLET", ha="center", fontsize=7.5, fontweight="bold", color=INK)

    def drawStats(fr):
        asx.clear()
        asx.axis("off")
        cur = [0.98]

        def line(txt, bold=False, mono=False, gap=0.0, color=MUTE):
            cur[0] -= gap
            asx.text(0.03, cur[0], txt, transform=asx.transAxes, va="top",
                      fontsize=(10 if bold else 8.6),
                      fontweight="bold" if bold else "normal",
                      family="monospace" if mono else "sans-serif",
                      color=INK if bold else color)
            cur[0] -= 0.046 * (txt.count("\n") + 1)

        line("Batch status", bold=True)
        cur[0] -= 0.01
        for st in sorted(fr["treeStatus"], key=lambda s: s["idx"]):
            i = st["idx"]
            ratio = ",".join(f"{v:.3f}" for v in ratioBySlot.get(i, []))
            status = "done" if st["done"] else ("loading" if st["active"] else "idle")
            line(f"slot{i+1}  d={depthBySlot.get(i,'?')}  {status}", gap=0.02)
            line(f"  target: {ratio}", mono=True, color=MUTE)

        cum = fr["cum"]
        line("Cumulative parallel K/B/L", bold=True, gap=0.035)
        line(f"K flows  : {cum[0]}\nB bends  : {cum[1]}\nL length : {cum[2]}",
             mono=True, gap=0.01)

        line("Reagent legend", bold=True, gap=0.035)
        seen = sorted({elm for row in fr["grid"] for elm in row
                        if isinstance(elm, str) and elm.startswith("R")},
                       key=lambda l: int(l[1:]))
        for lab in seen:
            asx.add_patch(mpatches.Circle((0.06, cur[0] - 0.012), 0.014,
                                            transform=asx.transAxes,
                                            facecolor=_reagentColor(lab),
                                            edgecolor="white", lw=0.6, zorder=3))
            asx.text(0.1, cur[0], lab, transform=asx.transAxes, va="top",
                      fontsize=8.4, color=MUTE)
            cur[0] -= 0.04

    def render(idx):
        fr = frames[idx]
        drawGrid(fr)
        drawStats(fr)

    if totCollisions:
        print(f"  ** animate_batch: {totCollisions} cross-tree cell collision(s) "
              f"blocked -- see [collision] log lines above **")

    anim = animation.FuncAnimation(fig, render, frames=len(frames), interval=1000 // fps)

    # ── GIF ──────────────────────────────────────────────────────────────────
    anim.save(outfile, writer=PillowWriter(fps=fps))
    print(f"  GIF saved : {outfile}")

    # ── MP4 video ─────────────────────────────────────────────────────────────
    # Uses imageio-ffmpeg (pip install imageio-ffmpeg) which bundles its own
    # ffmpeg binary -- no system PATH install required on Windows or Linux.
    import os as _os
    videoFile = _os.path.splitext(outfile)[0] + ".mp4"
    try:
        import numpy as _np
        import imageio_ffmpeg as _iff
        import imageio as _iio
        import io as _io
        _ffmpeg_exe = _iff.get_ffmpeg_exe()
        _frames_rgb = []
        for _idx in range(len(frames)):
            render(_idx)
            fig.canvas.draw()
            _buf = _io.BytesIO()
            fig.savefig(_buf, format='raw', dpi=fig.dpi)
            _buf.seek(0)
            _arr = _np.frombuffer(_buf.read(), dtype=_np.uint8)
            _w = int(fig.get_figwidth()  * fig.dpi)
            _h = int(fig.get_figheight() * fig.dpi)
            _frames_rgb.append(_arr.reshape(_h, _w, 4)[:, :, :3])
        _iio.mimwrite(videoFile, _frames_rgb, fps=fps,
                      ffmpeg_params=["-pix_fmt", "yuv420p",
                                     "-vcodec", "libx264",
                                     "-crf", "18"],
                      ffmpeg_log_level="quiet")
        print(f"  MP4 saved : {videoFile}")
    except ImportError:
        print("  MP4 skipped: run  pip install imageio[ffmpeg] imageio-ffmpeg")
    except Exception as _e:
        print(f"  MP4 skipped: {_e}")

    plt.close(fig)


# ----- 4. orchestrator -------------------------------------------------------
def prepareMultiTarget(targets, err, batch="batch", root="./output/",
                       baseIndex=0, render=True, animate=False, fps=10):
    """
    Generate a skewed tree for each target via repeated single-tree generation,
    place them together, and report combined preparation metrics.

    Every call is one *batch*: all of this batch's trees and results are written
    under  root/<batch>/  , each tree in its own  tr<N>/  subfolder.

    render=True writes a tree.png into each tree folder via graphviz; it needs
    the Graphviz system binaries (the `dot` command) on PATH. If `dot` is
    missing, rendering is skipped per tree without aborting the run.

    animate=True renders <batch>_sim.gif via animate_batch() (this file),
    requires matplotlib. Consistent by construction with the parallel K/B/L
    numbers, since both are built from the same buildTimelines() output.
    """
    batchDir = os.path.join(root, str(batch))
    os.makedirs(batchDir, exist_ok=True)
    if len(targets) > len(SLOTS):
        raise ValueError(
            f"At most {len(SLOTS)} targets fit on the {GRID_ROW}x{GRID_COL} grid"
        )

    # keep a copy of this batch's input targets alongside its trees
    with open(os.path.join(batchDir, "targets.csv"), "w", newline="") as tf:
        csv.writer(tf).writerows(targets)
    print(f"\n##### Batch '{batch}'  ->  {batchDir}  ({len(targets)} targets) #####")

    perTree = []
    allCells = []
    totWaste = totReagent = 0
    seqK = seqB = seqL = 0

    for i, ratio in enumerate(targets):
        ind = baseIndex + i
        treeDir = os.path.join(batchDir, f"tr{i + 1}")   # output/batch1/tr1, tr2, ...
        os.makedirs(treeDir, exist_ok=True)
        print(f"\n=== Target {i} -> {treeDir}: {ratio} (err={err}) ===")
        z3opFile, depth = generateSingleTree(ratio, err, ind, treeDir)
        if z3opFile is None:
            print(f"  no SAT skewed tree up to max depth; skipping target {i}")
            continue

        R = len(ratio)
        m = treeMetrics(z3opFile, depth, R, SLOTS[i])
        perTree.append((i, ratio, depth, m))

        # per-tree metrics, written inside that tree's own folder -- now
        # includes this tree's own BoundingBox/area, matching what
        # generateSkewedTree()'s single-target path gets from saveTree_getArea
        with open(os.path.join(treeDir, "metrics.csv"), "w", newline="") as mf:
            w = csv.writer(mf)
            w.writerow(["target", "depth", "waste", "reagent",
                        "BoundingBox", "area", "K", "B", "L"])
            w.writerow([",".join(map(str, ratio)), depth, m["waste"],
                        m["reagent"], m["boundingBox"], m["area"],
                        m["K"], m["B"], m["L"]])

        # render the tree image into this tree's folder (needs graphviz `dot`)
        if render:
            try:
                treeImg = os.path.join(treeDir, "tree")
                dot = cts.saveTree(m["x"], m["W"], treeImg)
                dot.render(treeImg, format="png", cleanup=True)
            except Exception as e:
                print(f"  (tree image skipped: {e})")

        allCells.extend(m["cells"])
        totWaste += m["waste"]
        totReagent += m["reagent"]
        seqK += m["K"]; seqB += m["B"]; seqL += m["L"]
        print(f"  depth={depth} waste={m['waste']} reagent={m['reagent']} "
              f"BoundingBox={m['boundingBox']} area={m['area']} "
              f"K={m['K']} B={m['B']} L={m['L']}")

    combinedArea = boundingBoxArea(allCells)
    footprint = len({tuple(c) for c in allCells})

    print("\n=== Combined multi-target preparation ===")
    print(f"trees placed      : {len(perTree)}")
    print(f"combined area     : {combinedArea}")
    print(f"footprint cells   : {footprint}")
    print(f"total waste       : {totWaste}")
    print(f"total reagent use : {totReagent}")
    print(f"sequential K/B/L  : {seqK} / {seqB} / {seqL}")

    par = parallelLoading(perTree)
    if par is not None:
        print(f"parallel   K/B/L  : {par[0]} / {par[1]} / {par[2]}")
        if par[3]:
            print(f"parallel   collisions blocked : {par[3]}")

    # log one summary row inside the batch folder
    resultFile = os.path.join(batchDir, f"multiResults_{err}.csv")
    with open(resultFile, "a", newline="") as op:
        w = csv.writer(op)
        row = [len(perTree), combinedArea, footprint, totWaste, totReagent,
               seqK, seqB, seqL]
        if par is not None:
            row += list(par)   # K, B, L, collisions
        w.writerow(row)

    # ---- graphical simulation of the whole batch on one chip --------------
    if animate and perTree:
        gifPath = os.path.join(batchDir, f"{batch}_sim.gif")
        animate_batch(perTree, outfile=gifPath, fps=fps)
        print(f"animation         : {gifPath}")

    return perTree


# ----- input helper ----------------------------------------------------------
def readTargets(path, drop_last=False):
    """
    Read target ratios, one per line, comma-separated floats.
    Set drop_last=True if your file follows the testCases_1.csv convention
    (ratios followed by one trailing column that should be discarded).
    """
    targets = []
    with open(path) as f:
        for line in csv.reader(f):
            vals = [float(v) for v in line if v != '']
            if not vals:
                continue
            targets.append(vals[:-1] if drop_last else vals)
    return targets


def runBatches(inputDir="./batches", err=0.001, root="./output/",
               drop_last=False, render=True, animate=False):
    """
    Treat every .csv in `inputDir` as a separate batch of targets. Each file
    'foo.csv' produces its own batch folder  root/foo/  containing that batch's
    tr<N>/ tree folders, targets.csv, and multiResults_<err>.csv.
    """
    if not os.path.isdir(inputDir):
        print(f"no batch directory '{inputDir}'")
        return
    files = sorted(f for f in os.listdir(inputDir) if f.lower().endswith(".csv"))
    if not files:
        print(f"no .csv batches found in '{inputDir}'")
        return
    for fname in files:
        batch = os.path.splitext(fname)[0]
        targets = readTargets(os.path.join(inputDir, fname), drop_last=drop_last)
        prepareMultiTarget(targets, err, batch=batch, root=root, render=render,
                           animate=animate)


def readRows(path, target_sep=';', drop_last=False):
    """
    Read a rows-input file. Each LINE is one batch (a "row"); within a line,
    targets are separated by `target_sep` (default ';') and the ratio components
    of one target are comma-separated.

        0.5,0.25,0.25;0.25,0.5,0.25;0.125,0.375,0.5   -> row with 3 targets
        0.625,0.125,0.25;0.375,0.375,0.25             -> row with 2 targets

    Returns a list of rows, where each row is a list of target ratios.
    Keep the reagent count consistent within a row (one chip reagent set).
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            targets = []
            for chunk in line.split(target_sep):
                chunk = chunk.strip()
                if not chunk:
                    continue
                vals = [float(v) for v in chunk.split(',') if v != '']
                if not vals:
                    continue
                targets.append(vals[:-1] if drop_last else vals)
            if targets:
                rows.append(targets)
    return rows


def runRows(path="input.csv", err=0.001, root="./output/", target_sep=';',
            drop_last=False, render=True, animate=False):
    """
    Run one batch per row of `path`. Row N -> root/rowN/tr1..trK .
    """
    rows = readRows(path, target_sep=target_sep, drop_last=drop_last)
    if not rows:
        print(f"no rows found in '{path}'")
        return
    print(f"{len(rows)} row(s) from '{path}'")
    for n, targets in enumerate(rows, start=1):
        prepareMultiTarget(targets, err, batch=f"row{n}", root=root,
                           render=render, animate=animate)


if __name__ == "__main__":
    err = 0.01

    if os.path.exists("input.csv"):
        # rows-input: each line is a batch -> output/row1, output/row2, ...
        runRows("input.csv", err=err, drop_last=False, animate=True)
    elif os.path.isdir("./batches"):
        # one batch folder per CSV file in ./batches/
        runBatches("./batches", err=err, drop_last=False, animate=True)
    elif os.path.exists("multiTargets.csv"):
        targets = readTargets("multiTargets.csv", drop_last=False)
        prepareMultiTarget(targets, err, batch="multiTargets", animate=True)
    else:
        # demo batch (requires z3 + NTM/LAFCADFL to actually run)
        targets = [
            [0.5,   0.25,  0.25],
            [0.25,  0.5,   0.25],
            [0.125, 0.375, 0.5],
        ]
        prepareMultiTarget(targets, err, batch="demo", animate=True)
