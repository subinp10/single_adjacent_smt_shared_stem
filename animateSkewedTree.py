"""
animateSkewedTree.py

Graphical simulation of ONE solved skewed tree using the REAL
LAFCADFL.getPlacementAndLoading pipeline -- no reimplementation.

The only modification to the real pipeline is:
  1. preStampGrid snapshot captured after LAFCA, before DFL wipes cells
  2. blockagePlacement globals read before they are reset, to identify
     which cells hold non-adjacent intermediates across timesteps

All cell assignments are done by LAFCA (createFile/Z3).
All routing is done by DFL.
findBlockages runs inside getPlacementAndLoading exactly as normal.

Modules required (must be on import path):
    createTreeForShared : getMixerData, getPlacement
    getLoadingData      : getLoadingData
    LAFCADFL            : getPlacementAndLoading, LAFCA_DFL globals,
                          LAFCA.createFile, DFL
matplotlib only imported if animateSkewedTree() is called.
"""

import createTreeForShared as cts
import getLoadingData as gld

# ── colours ───────────────────────────────────────────────────────────────────
_PAL   = ["#D85A30","#378ADD","#639922","#D4537E",
           "#EF9F27","#7F77DD","#1D9E75","#F09595"]
_MIXER = "#F2A623"
_TGT   = "#7F4FDB"
_EMPTY = "#f2f1ec"
_DOT   = "#d7dce4"
_INK   = "#2b2f36"

def _rc(label):
    try: idx = int(label[1:]) - 1
    except: idx = 0
    return _PAL[idx % len(_PAL)]

def _cc(elm):
    if not isinstance(elm, str) or elm == "*": return _EMPTY, ""
    if elm == "T":           return _TGT,   "T"
    if elm.startswith("R"): return _rc(elm), elm
    if elm.startswith("M"): return _MIXER,  elm
    return _EMPTY, ""


# ── timeline helpers (port of KBL() construction) ────────────────────────────
def buildParentAndTimestamps(M, assignment, depth):
    parent = {}
    for mix in M:
        for reagent in set(M[mix]):
            if not reagent.startswith('M'): continue
            i, j  = int(reagent[1:]), int(mix[1:])
            cnt   = M[mix].count(reagent)
            cells = assignment[mix][1]
            new_cells = ([cells[0], cells[2]] if cnt < 3 else
                         [cells[0], cells[1], cells[2]]) \
                        if i > j + 1 else cells[:]
            parent.setdefault(reagent, [])
            parent[reagent] += new_cells
    parent['M0'] = assignment['M0'][1][:]
    timeStampMap = {}
    for name, (ts, _) in assignment.items():
        timeStampMap.setdefault(ts, []).append(name)
    return parent, timeStampMap


def buildTimestepData(mixList, M, assignment):
    Mixtures = {}; loadingCells = {}
    reagentList = {}; blockageList = {}; units = {}
    for mix in mixList:
        cells = assignment[mix][1]
        Mixtures[mix]     = cells[:]
        loadingCells[mix] = cells[:]
        reagents = []; blockage = {}; unit = {}
        for reagent in M[mix]:
            if reagent.startswith('M'):
                if reagent not in blockage:
                    i, j = int(reagent[1:]), int(mix[1:])
                    cnt  = sum(1 for r in M[mix] if r == reagent)
                    unit[reagent] = cnt
                    if i > j + 1:
                        # Block left-column-first, leaving BR free -- matches
                        # the validated KBL() convention in LAFCA_DFL.py.
                        # cells ordering (from cts.getPlacement): TL,TR,BL,BR
                        # The non-adjacent block always sits to the left/above
                        # (getPlacement only ever shifts intervening mixers
                        # right or right+down), so keeping BR as the last
                        # free cell keeps a routable path open. Blocking
                        # TL+BL for cnt==2 is a vertically-adjacent pair
                        # (same column), not diagonal -- traceable for LAFCA.
                        if cnt == 1:   b = [cells[2]]                       # BL
                        elif cnt == 2: b = [cells[0], cells[2]]             # TL, BL
                        elif cnt == 3: b = [cells[0], cells[1], cells[2]]   # TL, TR, BL
                        blockage[reagent] = b
                    else:
                        blockage[reagent] = [c for c in assignment[reagent][1]
                                              if c in assignment[mix][1]]
            else:
                reagents.append(reagent)
        reagentList[mix]  = reagents
        blockageList[mix] = blockage
        units[mix]        = unit
    return Mixtures, loadingCells, reagentList, blockageList, units


# ── animation function ────────────────────────────────────────────────────────
def animateSkewedTree(z3opFile, R, d, startingCell, outfile,
                      fps=2, gridRow=10, gridCol=10):
    """
    Animate loading for one solved skewed tree.
    Uses the real LAFCADFL.getPlacementAndLoading directly.
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import matplotlib.animation as animation
    from matplotlib.animation import PillowWriter
    import os as _os, tempfile as _tmp

    # ── real LAFCADFL imports ─────────────────────────────────────────────────
    import LAFCADFL.LAFCA_DFL as _ld
    from LAFCADFL.LAFCA_DFL import getPlacementAndLoading
    from LAFCADFL.LAFCA import createFile

    # patch createFile to use absolute temp paths (Windows backslash fix)
    _tmpdir = _tmp.gettempdir().replace('\\', '/')
    _orig_createFile = createFile

    def _createFile_patched(mixtures, coordinates, infile, outfile_cf):
        _z3  = _os.path.join(_tmpdir, 'lafca_z3File.py').replace('\\', '/')
        _out = _os.path.join(_tmpdir, 'lafca_output.txt').replace('\\', '/')
        open(_out, 'w').close()
        return _orig_createFile(mixtures, coordinates, _z3, _out)

    import LAFCADFL.LAFCA as _lafca_mod
    _lafca_mod.createFile = _createFile_patched

    # ── tree data ─────────────────────────────────────────────────────────────
    Rvals, M, x, W = cts.getMixerData(z3opFile, R, d)
    assignment = cts.getPlacement(W, d, list(startingCell))
    K, B, L    = gld.getLoadingData(z3opFile, list(startingCell),
                                     d, gridRow, gridCol, 4, R)

    parent, timeStampMap = buildParentAndTimestamps(M, assignment, d)
    target = [v / (4 ** d) for v in Rvals[0]]
    inlet  = (0, gridRow - 1)
    outlet = (gridCol - 1, 0)

    # consumer timestep for each non-adjacent source mixer
    # {f'M{i}': consumerTs} -- used to set protection expiry
    consumerTs = {}
    for i in W:
        for j in W[i]:
            if i > j + 1:
                consumerTs[f'M{i}'] = d - 1 - j

    grid    = [['*'] * gridCol for _ in range(gridRow)]
    visGrid = [['*'] * gridCol for _ in range(gridRow)]
    frames  = []
    cumFlow = cumBend = cumLen = 0

    # protected: {(r,c): endTimestamp}
    # built live from findBlockages globals after each timestep
    protected = {}

    for t in sorted(timeStampMap):
        mixList = timeStampMap[t]
        Mixtures, loadingCells, reagentList, blockageList, units = \
            buildTimestepData(mixList, M, assignment)

        # ── snapshot AFTER LAFCA, BEFORE DFL ─────────────────────────────────
        # getPlacementAndLoading calls LAFCA then DFL -- DFL wipes the grid.
        # We need to see the grid between those two steps.
        # Strategy: wrap DFL to capture state just before it runs.
        _preStamp = [None]   # will hold the snapshot

        from LAFCADFL.DFL import DFL as _realDFL
        def _DFL_wrapper(inlet_arg, outlet_arg, g, reagents, cells):
            _preStamp[0] = [row[:] for row in g]   # snapshot here
            return _realDFL(inlet_arg, outlet_arg, g, reagents, cells)

        # patch DFL in the LAFCA_DFL namespace
        _ld.DFL = _DFL_wrapper

        # ── call the REAL getPlacementAndLoading ──────────────────────────────
        loadingPaths = getPlacementAndLoading(
            Mixtures, parent, loadingCells,
            reagentList, blockageList, units, grid
        )

        # restore real DFL
        _ld.DFL = _realDFL

        preStampGrid = _preStamp[0] if _preStamp[0] else [row[:] for row in grid]

        # ── read findBlockages result (globals already reset by now, but we
        #    can infer what was placed by comparing preStampGrid with grid)
        # cells that have an M-type label in preStampGrid but are NOT in
        # the current mixer's fresh reagents are blockage placements
        cellsThisStep = []
        for mix in mixList:
            cellsThisStep.extend(tuple(c) for c in assignment[mix][1])

        for (r, c) in cellsThisStep:
            val = preStampGrid[r][c]
            if isinstance(val, str) and val.startswith('M'):
                fluid = val
                endTs = consumerTs.get(fluid)
                if endTs is not None and endTs > t:
                    cur = protected.get((r, c))
                    if cur is None or endTs > cur:
                        protected[(r, c)] = endTs

        # ── cellLabel from LAFCA's actual output ──────────────────────────────
        cellLabel = {}
        for (r, c) in cellsThisStep:
            val = preStampGrid[r][c]
            visGrid[r][c] = '*'
            if isinstance(val, str) and val.startswith('R'):
                cellLabel[(r, c)] = val
            elif isinstance(val, str) and val.startswith('M'):
                visGrid[r][c] = val   # carried droplet: instant

        # ── stage 1: reveal flows progressively ──────────────────────────────
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
                               "paths": [(label, headPath)],
                               "t": t, "mixer": "+".join(mixList),
                               "cells": cellsThisStep,
                               "cum": (cumFlow, cumBend, cumLen)})
            for cell in path:
                tc = (cell[0], cell[1])
                if cellLabel.get(tc) == label:
                    visGrid[tc[0]][tc[1]] = label

        # ── stage 2: loaded ───────────────────────────────────────────────────
        frames.append({"grid": [row[:] for row in visGrid], "paths": [],
                       "t": t, "mixer": "+".join(mixList),
                       "cells": cellsThisStep,
                       "cum": (cumFlow, cumBend, cumLen)})

        # ── stage 3: mix-stamp (skip protected cells) ─────────────────────────
        for mix in mixList:
            stampLabel = "T" if int(mix[1:]) == 0 else mix
            for (r, c) in [tuple(c) for c in assignment[mix][1]]:
                endTs = protected.get((r, c))
                if endTs is not None and t < endTs:
                    continue   # non-adjacent intermediate must persist
                visGrid[r][c] = stampLabel
        frames.append({"grid": [row[:] for row in visGrid], "paths": [],
                       "t": t, "mixer": "+".join(mixList),
                       "cells": cellsThisStep,
                       "cum": (cumFlow, cumBend, cumLen)})

        # ── sync visGrid from real grid (skip protected cells) ────────────────
        for r in range(gridRow):
            for c in range(gridCol):
                endTs = protected.get((r, c))
                if endTs is not None and t < endTs:
                    continue
                visGrid[r][c] = grid[r][c]
        if "M0" in mixList:
            for (r, c) in [tuple(c) for c in assignment["M0"][1]]:
                visGrid[r][c] = "T"
        frames.append({"grid": [row[:] for row in visGrid], "paths": [],
                       "t": t, "mixer": "+".join(mixList),
                       "cells": cellsThisStep,
                       "cum": (cumFlow, cumBend, cumLen)})

    # restore patched functions
    _lafca_mod.createFile = _orig_createFile

    # ── render ────────────────────────────────────────────────────────────────
    def xy(r, c): return c, (gridRow - 1 - r)

    fig = plt.figure(figsize=(9, 6))
    fig.patch.set_facecolor("white")
    gs  = gridspec.GridSpec(1, 2, width_ratios=[2, 1], figure=fig, wspace=0.15)
    ag  = fig.add_subplot(gs[0, 0])
    asx = fig.add_subplot(gs[0, 1])

    def drawGrid(fr):
        ag.clear()
        ag.set_xlim(-0.7, gridCol - 0.3); ag.set_ylim(-0.7, gridRow - 0.3)
        ag.set_aspect("equal"); ag.axis("off")
        ag.text((gridCol-1)/2, gridRow-0.2,
                f"t = {fr['t']}  ({fr['mixer']})",
                ha="center", fontsize=11, fontweight="bold", color=_INK)
        for r in range(gridRow):
            for c in range(gridCol):
                px, py = xy(r, c)
                ag.add_patch(mpatches.Circle((px,py),0.07,
                    facecolor=_DOT,edgecolor="none",zorder=1))
        for r in range(gridRow):
            for c in range(gridCol):
                color, lbl = _cc(fr["grid"][r][c])
                if not lbl: continue
                px, py = xy(r, c)
                ag.add_patch(mpatches.FancyBboxPatch(
                    (px-0.42,py-0.42),0.84,0.84,
                    boxstyle="round,pad=0.02,rounding_size=0.12",
                    facecolor=color,edgecolor="white",linewidth=1.2,zorder=2))
                ag.text(px,py,lbl,ha="center",va="center",
                        fontsize=7,fontweight="bold",color="white",zorder=3)
        for lbl, path in fr["paths"]:
            if len(path) < 2: continue
            color = _rc(lbl)
            pts = [xy(p[0],p[1]) for p in path]
            ag.plot([p[0] for p in pts],[p[1] for p in pts],
                    color=color,lw=3.0,alpha=0.85,zorder=4,
                    solid_capstyle="round")
            hx,hy=pts[-1]
            ag.add_patch(mpatches.Circle((hx,hy),0.2,
                facecolor="#ffffff",edgecolor=color,lw=1.6,zorder=5))
        ix,iy=xy(*inlet); ox,oy=xy(*outlet)
        for (px,py),txt in [((ix,iy),"INLET"),((ox,oy),"OUTLET")]:
            ag.add_patch(mpatches.Circle((px,py),0.26,
                facecolor=_INK,edgecolor="none",zorder=5))
            dy = 0.5 if txt=="INLET" else -0.5
            ag.text(px,py+dy,txt,ha="center",fontsize=7,
                    fontweight="bold",color=_INK)
        rs=[c[0] for c in fr["cells"]]; cs=[c[1] for c in fr["cells"]]
        x0,y0=xy(max(rs),min(cs))
        ag.add_patch(plt.Rectangle((x0-0.46,y0-0.46),
            (max(cs)-min(cs))+0.92,(max(rs)-min(rs))+0.92,
            fill=False,edgecolor=_INK,lw=2.0,zorder=6))

    def drawStats(fr):
        asx.clear(); asx.axis("off")
        cur=[0.98]
        def line(txt,bold=False,mono=False,gap=0.0,color=None):
            cur[0]-=gap
            asx.text(0.03,cur[0],txt,transform=asx.transAxes,va="top",
                     fontsize=(10 if bold else 8.6),
                     fontweight="bold" if bold else "normal",
                     family="monospace" if mono else "sans-serif",
                     color=color or _INK)
            cur[0]-=0.05*(txt.count("\n")+1)
        line("Target ratio",bold=True)
        line(":".join(f"{v:.4f}" for v in target),mono=True,gap=0.02)
        line("K / B / L  (getLoadingData estimate)",bold=True,gap=0.04)
        line(f"K={K}  B={B}  L={L}",mono=True,gap=0.02)
        line("K / B / L  (live DFL)",bold=True,gap=0.04)
        cum=fr["cum"]
        line(f"K={cum[0]}  B={cum[1]}  L={cum[2]}",mono=True,gap=0.02)
        line("Active this timestep",bold=True,gap=0.04)
        line(f"depth={d}  mixer(s)={fr['mixer']}\ncells={fr['cells']}",
             mono=True,gap=0.02)

    def render(idx):
        fr=frames[idx]; drawGrid(fr); drawStats(fr)

    anim=animation.FuncAnimation(fig,render,frames=len(frames),
                                  interval=1000//fps)
    anim.save(outfile,writer=PillowWriter(fps=fps))
    print(f"  GIF saved : {outfile}")

    videoFile=_os.path.splitext(outfile)[0]+".mp4"
    try:
        import numpy as _np, imageio_ffmpeg as _iff
        import imageio as _iio, io as _io
        _iff.get_ffmpeg_exe()
        _frames=[]
        for _i in range(len(frames)):
            render(_i); fig.canvas.draw()
            _buf=_io.BytesIO(); fig.savefig(_buf,format='raw',dpi=fig.dpi)
            _buf.seek(0)
            _arr=_np.frombuffer(_buf.read(),dtype=_np.uint8)
            _w=int(fig.get_figwidth()*fig.dpi); _h=int(fig.get_figheight()*fig.dpi)
            _frames.append(_arr.reshape(_h,_w,4)[:,:,:3])
        _iio.mimwrite(videoFile,_frames,fps=fps,
                      ffmpeg_params=["-pix_fmt","yuv420p","-vcodec","libx264","-crf","18"],
                      ffmpeg_log_level="quiet")
        print(f"  MP4 saved : {videoFile}")
    except ImportError:
        print("  MP4 skipped: pip install imageio[ffmpeg] imageio-ffmpeg")
    except Exception as _e:
        print(f"  MP4 skipped: {_e}")
    plt.close(fig)
