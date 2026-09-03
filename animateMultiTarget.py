"""
animateMultiTarget.py

Graphical simulation for multiTargetLAFCA.py's TRUE parallel placement +
routing pipeline. Mirrors how animateSkewedTree.py sits alongside stg1.py
in the single-target project: the core pipeline (multiTargetLAFCA.py)
imports and calls this module's animateMultiTarget() when animate=True,
rather than the animation code living inline.

Uses the SAME real pipeline as multiTargetLAFCA()'s non-animated path:
    placement : buildTimestepData + placeOneTreeTimestep, called ONCE PER
                TREE per timestep (LAFCA Z3 cell assignment + the
                corrected findBlockages/dfs corner selection)
    routing   : DFL, called EXACTLY ONCE PER TIMESTEP across ALL active
                trees' combined cells-to-load -- true parallel routing,
                not one DFL call per tree.

Format matches animateSkewedTree_corrected.py: 2-panel layout (grid +
stats sidebar), progressive reagent-label reveal tied to path progress
(reagents appear only as their route reaches them, not immediately on
placement), and the root (M0) shown as "T" in the display only -- the
real grid keeps "M0" since other logic depends on that bookkeeping.
Per-tree colored box outlines distinguish which tree's mixer is which,
since (unlike the single-target case) several trees are active at once.
"""

from LAFCADFL.DFL import DFL as _realDFL
from multiTargetLAFCA import (
    buildTimestepData,
    placeOneTreeTimestep,
    relabelAfterRouting,
    GRID_ROW,
    GRID_COL,
)

# ----- colors ----------------------------------------------------------------
_PAL = ["#D85A30", "#378ADD", "#639922", "#D4537E",
        "#EF9F27", "#7F77DD", "#1D9E75", "#F09595"]
_MIXER = "#F2A623"
_TARGET = "#7F4FDB"
_EMPTY = "#f2f1ec"
_DOT = "#d7dce4"
_INK = "#2b2f36"
_BOX_COLORS = ["#2b2f36", "#993c1d", "#3c3489", "#085041",
               "#712b13", "#26215c", "#04342c", "#4a1b0c", "#4b1528"]


def _reagentColor(label):
    try:
        idx = int(label[1:]) - 1
    except (ValueError, IndexError):
        idx = 0
    return _PAL[idx % len(_PAL)]


def _cellColor(elm):
    if not isinstance(elm, str) or elm == "*":
        return _EMPTY, ""
    if elm == "T":
        return _TARGET, "T"
    if elm.startswith("R"):
        return _reagentColor(elm), elm
    if elm.startswith("M"):
        return _MIXER, elm
    return _EMPTY, ""


def animateMultiTarget(perTree, outfile, fps=2,
                        grid_row=GRID_ROW, grid_col=GRID_COL):
    """
    Runs the real placement + parallel-routing pipeline for every tree in
    perTree (as built by multiTargetLAFCA()), capturing a frame after
    each timestep's placement, after each routed flow (progressively
    revealed), and after each timestep's post-routing relabel/wash --
    then renders and saves a GIF.

    perTree: list of dicts, one per tree, each with keys
        "ratio", "depth", "M", "assignment", "parent", "timeStampMap",
        and optionally "estK"/"estB"/"estL" (getLoadingData estimate,
        shown in the sidebar alongside the live DFL numbers).

    Returns (totalFlow, totalBendings, totalPathLength).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import matplotlib.animation as animation
    from matplotlib.animation import PillowWriter

    inlet, outlet = [0, grid_row - 1], [grid_col - 1, 0]
    grid = [['*'] * grid_col for _ in range(grid_row)]      # real: keeps 'M0' etc,
                                                              # fed to placement/DFL
    visGrid = [['*'] * grid_col for _ in range(grid_row)]   # display: 'M0' -> 'T'
    frames = []
    cumFlow = cumBend = cumLen = 0
    maxT = max(d["depth"] for d in perTree)

    estK = sum(d.get("estK", 0) for d in perTree)
    estB = sum(d.get("estB", 0) for d in perTree)
    estL = sum(d.get("estL", 0) for d in perTree)

    for t in range(maxT):
        combinedCellsToLoad = []
        combinedReagents = set()
        treeMixturesThisT = []
        activeBoxesThisT = []   # (boxCells, colorIdx, treeIdx, mixName)
        cellLabel = {}          # (r,c) -> expected fresh-reagent label, for reveal

        for treeIdx, d in enumerate(perTree):
            if t not in d["timeStampMap"]:
                continue
            mixList = d["timeStampMap"][t]
            Mixtures, loadingCells, reagentList, blockageList, units = \
                buildTimestepData(mixList, d["M"], d["assignment"])

            cellsToLoad, reagents = placeOneTreeTimestep(
                loadingCells, reagentList, blockageList, units, grid)
            combinedCellsToLoad.extend(cellsToLoad)
            combinedReagents.update(reagents)
            treeMixturesThisT.append((Mixtures, d["parent"]))

            colorIdx = treeIdx % len(_BOX_COLORS)
            for mix in mixList:
                boxCells = [tuple(c) for c in d["assignment"][mix][1]]
                activeBoxesThisT.append((boxCells, colorIdx, treeIdx, mix))
                for (r, c) in boxCells:
                    val = grid[r][c]
                    visGrid[r][c] = '*'
                    if isinstance(val, str) and val.startswith('R'):
                        cellLabel[(r, c)] = val
                    elif isinstance(val, str) and val.startswith('M'):
                        visGrid[r][c] = val   # carried droplet: instant, not routed

        frames.append({"grid": [row[:] for row in visGrid], "paths": [], "t": t,
                        "boxes": activeBoxesThisT,
                        "cum": (cumFlow, cumBend, cumLen), "est": (estK, estB, estL)})

        # ---- ONE DFL call for this timestep, across ALL active trees ----
        if combinedCellsToLoad:
            loadingPaths = _realDFL(inlet, outlet, grid,
                                     list(combinedReagents), combinedCellsToLoad)
            for order in loadingPaths:
                label, bends, path = order[0], order[1], order[2]
                cumFlow += 1; cumBend += bends; cumLen += len(path)
                step = max(1, len(path) // 3)
                for kk in range(step, len(path) + 1, step):
                    headPath = path[:kk]
                    for cell in headPath:
                        tc = (cell[0], cell[1])
                        if cellLabel.get(tc) == label:
                            visGrid[tc[0]][tc[1]] = label
                    frames.append({"grid": [row[:] for row in visGrid],
                                   "paths": [(label, headPath)], "t": t,
                                   "boxes": activeBoxesThisT,
                                   "cum": (cumFlow, cumBend, cumLen),
                                   "est": (estK, estB, estL)})
                for cell in path:
                    tc = (cell[0], cell[1])
                    if cellLabel.get(tc) == label:
                        visGrid[tc[0]][tc[1]] = label
                frames.append({"grid": [row[:] for row in visGrid],
                               "paths": [(label, path)], "t": t,
                               "boxes": activeBoxesThisT,
                               "cum": (cumFlow, cumBend, cumLen),
                               "est": (estK, estB, estL)})

        # relabel/wash every active tree's box AFTER routing -- real grid
        # keeps 'M0', visGrid shows 'T' for the root
        for Mixtures, parentMix in treeMixturesThisT:
            relabelAfterRouting(Mixtures, parentMix, grid, visGrid)

        frames.append({"grid": [row[:] for row in visGrid], "paths": [], "t": t,
                        "boxes": activeBoxesThisT,
                        "cum": (cumFlow, cumBend, cumLen), "est": (estK, estB, estL)})

    # ----- render: 2-panel, grid + stats sidebar -----
    def xy(r, c):
        return c, (grid_row - 1 - r)

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
                 ha="center", fontsize=11, fontweight="bold", color=_INK)

        for r_ in range(grid_row):
            for c_ in range(grid_col):
                px, py = xy(r_, c_)
                ag.add_patch(mpatches.Circle((px, py), 0.07,
                                              facecolor=_DOT, edgecolor="none", zorder=1))

        for r_ in range(grid_row):
            for c_ in range(grid_col):
                color, lbl = _cellColor(fr["grid"][r_][c_])
                if not lbl:
                    continue
                px, py = xy(r_, c_)
                ag.add_patch(mpatches.FancyBboxPatch(
                    (px - 0.42, py - 0.42), 0.84, 0.84,
                    boxstyle="round,pad=0.02,rounding_size=0.12",
                    facecolor=color, edgecolor="white", linewidth=1.2, zorder=2))
                ag.text(px, py, lbl, ha="center", va="center",
                        fontsize=7, fontweight="bold", color="white", zorder=3)

        for lbl, path in fr["paths"]:
            if len(path) < 2:
                continue
            color = _reagentColor(lbl)
            pts = [xy(p[0], p[1]) for p in path]
            ag.plot([p[0] for p in pts], [p[1] for p in pts],
                    color=color, lw=3.0, alpha=0.85, zorder=4, solid_capstyle="round")

        # per-tree colored box outline
        for boxCells, colorIdx, treeIdx, mixName in fr["boxes"]:
            rs = [c[0] for c in boxCells]
            cs = [c[1] for c in boxCells]
            x0, y0 = xy(max(rs), min(cs))
            ag.add_patch(plt.Rectangle((x0 - 0.46, y0 - 0.46),
                         (max(cs) - min(cs)) + 0.92, (max(rs) - min(rs)) + 0.92,
                         fill=False, edgecolor=_BOX_COLORS[colorIdx], lw=2.2, zorder=6))

        ix, iy = xy(*inlet)
        ox, oy = xy(*outlet)
        for (px, py), txt in [((ix, iy), "INLET"), ((ox, oy), "OUTLET")]:
            ag.add_patch(mpatches.Circle((px, py), 0.26, facecolor=_INK,
                                          edgecolor="none", zorder=7))
            dy = 0.5 if txt == "INLET" else -0.5
            ag.text(px, py + dy, txt, ha="center", fontsize=7,
                    fontweight="bold", color=_INK)

    def drawStats(fr):
        asx.clear(); asx.axis("off")
        cur = [0.98]

        def line(txt, bold=False, mono=False, gap=0.0, color=None):
            cur[0] -= gap
            asx.text(0.03, cur[0], txt, transform=asx.transAxes, va="top",
                      fontsize=(10 if bold else 8.6),
                      fontweight="bold" if bold else "normal",
                      family="monospace" if mono else "sans-serif",
                      color=color or _INK)
            cur[0] -= 0.045 * (txt.count("\n") + 1)

        line("Targets", bold=True)
        for i, d in enumerate(perTree):
            ratioStr = ":".join(f"{v:.3f}".rstrip('0').rstrip('.') or "0"
                                 for v in d["ratio"])
            color = _BOX_COLORS[i % len(_BOX_COLORS)]
            line(f"tr{i + 1}: {ratioStr}  (depth={d['depth']})",
                 mono=True, gap=0.012, color=color)

        line("K / B / L  (getLoadingData estimate)", bold=True, gap=0.04)
        est = fr["est"]
        line(f"K={est[0]}  B={est[1]}  L={est[2]}", mono=True, gap=0.02)

        line("K / B / L  (live, parallel DFL)", bold=True, gap=0.04)
        cum = fr["cum"]
        line(f"K={cum[0]}  B={cum[1]}  L={cum[2]}", mono=True, gap=0.02)

        line("Active this timestep", bold=True, gap=0.04)
        active = [f"tr{ti + 1}:{mn}" for (_, _, ti, mn) in fr["boxes"]]
        line(f"t={fr['t']}\n{', '.join(active) if active else '(none)'}",
             mono=True, gap=0.02)

    def render(idx):
        fr = frames[idx]
        drawGrid(fr)
        drawStats(fr)

    anim = animation.FuncAnimation(fig, render, frames=len(frames),
                                    interval=1000 // fps)
    anim.save(outfile, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(f"  GIF saved : {outfile}")
    print(f"Parallel Flow: {cumFlow}, Parallel Bendings: {cumBend}, "
          f"Parallel Path Length: {cumLen}")
    return cumFlow, cumBend, cumLen
