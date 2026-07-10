"""
createTreeSharedStem.py

Draws the shared-stem tree from a result dict returned by
singleSMT.generateSharedStemTreeSingleSMT().

Since the combined approach produces ONE z3 output file (not two),
this reads stem + T1 + T2 all from that single file.

Structure:
    Stem   : M1 (leaf) ... Mc (root/split)
    T1     : A_(c+1) (leaf, gets w1) ... A_(d1) (root = Target1)
    T2     : B_(c+1) (leaf, gets w2) ... B_(d2) (root = Target2)

Shows mixer nodes with their mixing ratio, reagent nodes, and edge weights.
"""

import os


def _parseZ3(z3File):
    """Parse name=value pairs from a z3 output file into a dict of ints."""
    vals = {}
    skip = ('t_', 'At_', 'Bt_', 'tT_', 'tU_', 'S1', 'S2')
    with open(z3File, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if ' = ' not in line:
                continue
            name, val = line.split(' = ', 1)
            name = name.strip()
            if name.startswith(skip):
                continue
            try:
                vals[name] = int(val.strip())
            except ValueError:
                pass
    return vals


def drawSharedStemTree(result, outPath, R):
    """
    Draw the tree from a generateSharedStemTreeSingleSMT() result dict.

    Parameters
    ----------
    result  : dict returned by generateSharedStemTreeSingleSMT
              must contain 'c', 'd1', 'd2', 'z3file'
    outPath : output path without .png extension
    R       : number of reagents
    """
    try:
        import graphviz
    except ImportError:
        print("pip install graphviz")
        return None

    c  = result['c']
    d1 = result['d1']
    d2 = result['d2']

    v = _parseZ3(result['z3file'])   # stem + T1 + T2, all in one file

    dot = graphviz.Digraph(comment="SharedStemTree")
    dot.attr(rankdir='BT')
    dot.attr('graph', fontname='Helvetica', fontsize='11')
    dot.attr('node',  fontname='Helvetica', fontsize='11')
    dot.attr('edge',  fontname='Helvetica', fontsize='10')

    # -- stem nodes M1..Mc --
    for i in range(1, c + 1):
        color = 'orange' if i == c else 'lightyellow'
        Rs = [v.get(f"R_{i}_{r}", 0) for r in range(1, R + 1)]
        denom = 4 ** i
        ratio_str = ":".join(f"{val/denom:.3f}" for val in Rs)
        label = f"M{i}\n{ratio_str}"
        dot.node(f"M{i}", label, shape='ellipse',
                 style='filled', fillcolor=color)

    for i in range(2, c + 1):
        W = v.get(f"W_{i}_{i-1}", "?")
        dot.edge(f"M{i-1}", f"M{i}", label=str(W))

    for i in range(1, c + 1):
        for r in range(1, R + 1):
            val = v.get(f"x_{i}_{r}", 0)
            if val > 0:
                nid = f"sR{i}_{r}"
                dot.node(nid, f"R{r}", shape='ellipse',
                         style='filled', fillcolor='white')
                dot.edge(nid, f"M{i}", label=str(val))

    # -- T1 chain nodes A_(c+1)..A_(d1) --
    w1 = v.get('w1', '?')
    for i in range(c + 1, d1 + 1):
        As = [v.get(f"A_{i}_{r}", 0) for r in range(1, R + 1)]
        denom = 4 ** i
        ratio_str = ":".join(f"{val/denom:.3f}" for val in As)
        tag = "  (T1 target)" if i == d1 else ""
        label = f"A{i}\n{ratio_str}{tag}"
        dot.node(f"A{i}", label, shape='ellipse',
                 style='filled', fillcolor='lightblue')

    dot.edge(f"M{c}", f"A{c+1}", label=str(w1),
             style='dashed', color='steelblue')

    for i in range(c + 2, d1 + 1):
        AW = v.get(f"AW_{i}_{i-1}", "?")
        dot.edge(f"A{i-1}", f"A{i}", label=str(AW))

    for i in range(c + 1, d1 + 1):
        for r in range(1, R + 1):
            val = v.get(f"Ax_{i}_{r}", 0)
            if val > 0:
                nid = f"tR{i}_{r}"
                dot.node(nid, f"R{r}", shape='ellipse',
                         style='filled', fillcolor='white')
                dot.edge(nid, f"A{i}", label=str(val))

    # -- T2 chain nodes B_(c+1)..B_(d2) --
    w2 = v.get('w2', '?')
    for i in range(c + 1, d2 + 1):
        Bs = [v.get(f"B_{i}_{r}", 0) for r in range(1, R + 1)]
        denom = 4 ** i
        ratio_str = ":".join(f"{val/denom:.3f}" for val in Bs)
        tag = "  (T2 target)" if i == d2 else ""
        label = f"B{i}\n{ratio_str}{tag}"
        dot.node(f"B{i}", label, shape='ellipse',
                 style='filled', fillcolor='lightgreen')

    dot.edge(f"M{c}", f"B{c+1}", label=str(w2),
             style='dashed', color='darkgreen')

    for i in range(c + 2, d2 + 1):
        BW = v.get(f"BW_{i}_{i-1}", "?")
        dot.edge(f"B{i-1}", f"B{i}", label=str(BW))

    for i in range(c + 1, d2 + 1):
        for r in range(1, R + 1):
            val = v.get(f"Bx_{i}_{r}", 0)
            if val > 0:
                nid = f"uR{i}_{r}"
                dot.node(nid, f"R{r}", shape='ellipse',
                         style='filled', fillcolor='white')
                dot.edge(nid, f"B{i}", label=str(val))

    outPath = outPath.replace('\\', '/')
    dot.render(outPath, format='png', cleanup=True)
    print(f"Tree saved: {outPath}.png")
    return dot


if __name__ == "__main__":
    from singleSMT import generateSharedStemTreeSingleSMT

    result = generateSharedStemTreeSingleSMT(
        target1=[0.25, 0.50, 0.25],
        target2=[0.50, 0.25, 0.25],
        err=0.01,
        outDir="./sharedStem/",
        label="draw",
        maxDepth=13,
    )

    if result and "note" not in result:
        drawSharedStemTree(result, "./sharedStem/sharedTree", R=3)
