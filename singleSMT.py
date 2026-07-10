"""
singleSMT.py

Shared-stem skewed mixing tree for TWO target ratios -- COMBINED SMT
(Idea 1: stem + T1 chain + T2 chain solved TOGETHER in one Z3 instance).

Structure:
    Stem   : i=1 (leaf, pure fresh) -> i=c (root = split point)
    T1     : i=c+1 (leaf, gets w1 from Mc) -> i=d1 (root = Target1)
    T2     : i=c+1 (leaf, gets w2 from Mc) -> i=d2 (root = Target2)

Unlike the enumerate-then-test approach (Option B), here R_c_j is a
Z3 VARIABLE (not a fixed constant) referenced by BOTH the T1 leaf
equation and the T2 leaf equation. Sharing is guaranteed by construction:
there is only one R_c_j per reagent, used by both chains simultaneously.

Search:
    for c = c_start downto 1:
        for d1 = c+1 .. maxDepth:
            for d2 = c+1 .. maxDepth:
                build ONE combined SMT (stem + T1 + T2)
                solve
                if SAT: return result
"""

import os, sys, subprocess

N          = 4   # 2x2 mixer = 4 cells
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# UTILITY HELPERS
# =============================================================================

def _w(opfile, line=""):
    opfile.write(line + "\n")


def _declInts(opfile, *names):
    """Int() for a single variable, Ints() for multiple."""
    if len(names) == 1:
        _w(opfile, f"{names[0]} = Int('{names[0]}')")
    else:
        _w(opfile, ", ".join(names) + " = Ints('" + " ".join(names) + "')")


def _isSAT(path):
    if not os.path.exists(path):
        return False
    with open(path, encoding='utf-8') as f:
        return " = " in f.read()


def _parseZ3(path):
    vals = {}
    with open(path, encoding='utf-8') as f:
        for line in f:
            if " = " in line:
                k, v = line.strip().split(" = ", 1)
                try:
                    vals[k.strip()] = int(v.strip())
                except ValueError:
                    pass
    return vals


def _findDepth(target, err, R, maxDepth=10):
    """Minimum depth d such that target is representable at denom 4^d."""
    for d in range(1, maxDepth + 1):
        denom = N ** d
        vals  = [target[r] * denom for r in range(R)]
        if all(abs(v - round(v)) < 1e-9 for v in vals):
            rounded = [round(v) for v in vals]
            if all(abs(rounded[r] / denom - target[r]) <= err for r in range(R)):
                return d
    return maxDepth


# =============================================================================
# MODULE 1 -- SCRIPT HEADER
# =============================================================================

def initOPT(opfile):
    _w(opfile, "import sys, time")
    _w(opfile, "from z3 import *")
    _w(opfile, "")
    _w(opfile, "s = Optimize()")
    _w(opfile, "")


def finishOPT(opfile, z3opFile):
    z3opFile = z3opFile.replace("\\", "/")
    _w(opfile, "")
    _w(opfile, "s.set('timeout', 7200000)")
    _w(opfile, "t0     = time.time()")
    _w(opfile, "result = s.check()")
    _w(opfile, f"fp = open('{z3opFile}', 'w', encoding='utf-8')")
    _w(opfile, "if result == sat:")
    _w(opfile, "    m = s.model()")
    _w(opfile, "    for d in m:")
    _w(opfile, "        fp.write(str(d) + ' = ' + str(m[d]) + '\\n')")
    _w(opfile, "else:")
    _w(opfile, "    fp.write('unsat')")
    _w(opfile, "fp.flush()")
    _w(opfile, "fp.close()")
    _w(opfile, "print('Time =', round(time.time()-t0, 2), 's')")


# =============================================================================
# MODULE 2 -- VARIABLE DECLARATIONS (stem + T1 + T2, all in one script)
# =============================================================================

def declareVariables(c, d1, d2, R, opfile):
    """
    Declare every variable needed for the combined model.

    Stem (i=1..c):
        R_i_j, x_i_j, y_i_j, W_i_(i-1), t_i_(i-1)_r

    Split (shared by both chains, references the SAME R_c):
        w1, w2
        tT_r = w1 * R_c_r     (T1's use of the shared stem output)
        tU_r = w2 * R_c_r     (T2's use of the shared stem output)

    T1 chain (i=c+1..d1):
        A_i_j, Ax_i_j, Ay_i_j, AW_i_(i-1), At_i_(i-1)_r

    T2 chain (i=c+1..d2):
        B_i_j, Bx_i_j, By_i_j, BW_i_(i-1), Bt_i_(i-1)_r
    """
    _w(opfile, "# --- Stem variables (i=1..c) ---")
    for i in range(1, c + 1):
        _declInts(opfile, *[f"R_{i}_{j}" for j in range(1, R + 1)])
    for i in range(1, c + 1):
        _declInts(opfile, *[f"x_{i}_{j}" for j in range(1, R + 1)])
    for i in range(1, c + 1):
        _declInts(opfile, *[f"y_{i}_{j}" for j in range(1, R + 1)])
    for i in range(2, c + 1):
        _declInts(opfile, f"W_{i}_{i-1}")
    for i in range(2, c + 1):
        _declInts(opfile, *[f"t_{i}_{i-1}_{r}" for r in range(1, R + 1)])

    _w(opfile, "# --- Split (w1, w2 both reference the SAME R_c) ---")
    _declInts(opfile, "w1", "w2")
    _declInts(opfile, *[f"tT_{r}" for r in range(1, R + 1)])
    _declInts(opfile, *[f"tU_{r}" for r in range(1, R + 1)])

    _w(opfile, "# --- T1 chain variables (i=c+1..d1) ---")
    for i in range(c + 1, d1 + 1):
        _declInts(opfile, *[f"A_{i}_{j}" for j in range(1, R + 1)])
    for i in range(c + 1, d1 + 1):
        _declInts(opfile, *[f"Ax_{i}_{j}" for j in range(1, R + 1)])
    for i in range(c + 1, d1 + 1):
        _declInts(opfile, *[f"Ay_{i}_{j}" for j in range(1, R + 1)])
    for i in range(c + 2, d1 + 1):
        _declInts(opfile, f"AW_{i}_{i-1}")
    for i in range(c + 2, d1 + 1):
        _declInts(opfile, *[f"At_{i}_{i-1}_{r}" for r in range(1, R + 1)])

    _w(opfile, "# --- T2 chain variables (i=c+1..d2) ---")
    for i in range(c + 1, d2 + 1):
        _declInts(opfile, *[f"B_{i}_{j}" for j in range(1, R + 1)])
    for i in range(c + 1, d2 + 1):
        _declInts(opfile, *[f"Bx_{i}_{j}" for j in range(1, R + 1)])
    for i in range(c + 1, d2 + 1):
        _declInts(opfile, *[f"By_{i}_{j}" for j in range(1, R + 1)])
    for i in range(c + 2, d2 + 1):
        _declInts(opfile, f"BW_{i}_{i-1}")
    for i in range(c + 2, d2 + 1):
        _declInts(opfile, *[f"Bt_{i}_{i-1}_{r}" for r in range(1, R + 1)])
    _w(opfile, "")


# =============================================================================
# MODULE 3 -- LINEARITY CONSTRAINTS
# =============================================================================

def linearityConstraints(c, d1, d2, R, opfile):
    """
    Eliminate all products of two unknown variables via case split.

    Stem:  t_i_(i-1)_r = W_i_(i-1) * R_(i-1)_r        (W in {1,2,3})
    Split: tT_r = w1 * R_c_r,  tU_r = w2 * R_c_r      (w1,w2 in {1,2})
           NOTE: R_c_r is now a VARIABLE (not a constant like in Option B),
           so this product genuinely needs the case-split trick here.
    T1:    At_i_(i-1)_r = AW_i_(i-1) * A_(i-1)_r      (AW in {1,2,3})
    T2:    Bt_i_(i-1)_r = BW_i_(i-1) * B_(i-1)_r      (BW in {1,2,3})
    """
    _w(opfile, "# --- Stem linearity ---")
    for i in range(2, c + 1):
        for r in range(1, R + 1):
            for v in [1, 2, 3]:
                _w(opfile,
                   f"s.add(Implies(W_{i}_{i-1}=={v},"
                   f" t_{i}_{i-1}_{r}=={v}*R_{i-1}_{r}))")

    _w(opfile, "# --- Split linearity (R_c is a VARIABLE here) ---")
    for r in range(1, R + 1):
        for v in [1, 2]:
            _w(opfile, f"s.add(Implies(w1=={v}, tT_{r}=={v}*R_{c}_{r}))")
            _w(opfile, f"s.add(Implies(w2=={v}, tU_{r}=={v}*R_{c}_{r}))")

    _w(opfile, "# --- T1 chain linearity ---")
    for i in range(c + 2, d1 + 1):
        for r in range(1, R + 1):
            for v in [1, 2, 3]:
                _w(opfile,
                   f"s.add(Implies(AW_{i}_{i-1}=={v},"
                   f" At_{i}_{i-1}_{r}=={v}*A_{i-1}_{r}))")

    _w(opfile, "# --- T2 chain linearity ---")
    for i in range(c + 2, d2 + 1):
        for r in range(1, R + 1):
            for v in [1, 2, 3]:
                _w(opfile,
                   f"s.add(Implies(BW_{i}_{i-1}=={v},"
                   f" Bt_{i}_{i-1}_{r}=={v}*B_{i-1}_{r}))")
    _w(opfile, "")


# =============================================================================
# MODULE 4 -- DOMAIN AND CELL BUDGET CONSTRAINTS
# =============================================================================

def domainAndBudget(c, d1, d2, R, opfile):
    """
    Domain bounds:  x/Ax/Bx in [0,3],  W/AW/BW in [1,3],  w1/w2 in [1,2]
    Cell budgets:   every mixer's 4 cells = fresh + incoming intermediate
    """
    _w(opfile, "# --- Stem domain ---")
    for i in range(1, c + 1):
        for j in range(1, R + 1):
            _w(opfile, f"s.add(x_{i}_{j}>=0, x_{i}_{j}<=3)")
            _w(opfile, f"s.add(R_{i}_{j}>=0)")
    for i in range(2, c + 1):
        _w(opfile, f"s.add(W_{i}_{i-1}>=1, W_{i}_{i-1}<=3)")

    _w(opfile, "# --- Split domain ---")
    _w(opfile, "s.add(w1>=1, w1<=2)")
    _w(opfile, "s.add(w2>=1, w2<=2)")

    _w(opfile, "# --- T1 chain domain ---")
    for i in range(c + 1, d1 + 1):
        for j in range(1, R + 1):
            _w(opfile, f"s.add(Ax_{i}_{j}>=0, Ax_{i}_{j}<=3)")
            _w(opfile, f"s.add(A_{i}_{j}>=0)")
    for i in range(c + 2, d1 + 1):
        _w(opfile, f"s.add(AW_{i}_{i-1}>=1, AW_{i}_{i-1}<=3)")

    _w(opfile, "# --- T2 chain domain ---")
    for i in range(c + 1, d2 + 1):
        for j in range(1, R + 1):
            _w(opfile, f"s.add(Bx_{i}_{j}>=0, Bx_{i}_{j}<=3)")
            _w(opfile, f"s.add(B_{i}_{j}>=0)")
    for i in range(c + 2, d2 + 1):
        _w(opfile, f"s.add(BW_{i}_{i-1}>=1, BW_{i}_{i-1}<=3)")

    _w(opfile, "# --- Stem cell budgets ---")
    _w(opfile, "# Leaf i=1: sum(x_1_j) = 4 (pure fresh)")
    terms = "+".join(f"x_1_{r}" for r in range(1, R + 1))
    _w(opfile, f"s.add({terms}=={N})")
    _w(opfile, "# Inner i=2..c: sum(x_i_j) + W_(i,i-1) = 4")
    for i in range(2, c + 1):
        terms = "+".join(f"x_{i}_{r}" for r in range(1, R + 1))
        _w(opfile, f"s.add({terms}+W_{i}_{i-1}=={N})")

    _w(opfile, "# --- T1 cell budgets ---")
    _w(opfile, "# Leaf i=c+1: sum(Ax_(c+1)_j) + w1 = 4")
    terms = "+".join(f"Ax_{c+1}_{r}" for r in range(1, R + 1))
    _w(opfile, f"s.add({terms}+w1=={N})")
    _w(opfile, "# Inner i=c+2..d1: sum(Ax_i_j) + AW_(i,i-1) = 4")
    for i in range(c + 2, d1 + 1):
        terms = "+".join(f"Ax_{i}_{r}" for r in range(1, R + 1))
        _w(opfile, f"s.add({terms}+AW_{i}_{i-1}=={N})")

    _w(opfile, "# --- T2 cell budgets ---")
    _w(opfile, "# Leaf i=c+1: sum(Bx_(c+1)_j) + w2 = 4")
    terms = "+".join(f"Bx_{c+1}_{r}" for r in range(1, R + 1))
    _w(opfile, f"s.add({terms}+w2=={N})")
    _w(opfile, "# Inner i=c+2..d2: sum(Bx_i_j) + BW_(i,i-1) = 4")
    for i in range(c + 2, d2 + 1):
        terms = "+".join(f"Bx_{i}_{r}" for r in range(1, R + 1))
        _w(opfile, f"s.add({terms}+BW_{i}_{i-1}=={N})")

    _w(opfile, "# --- y flags (used in minimize) ---")
    for i in range(1, c + 1):
        for r in range(1, R + 1):
            _w(opfile, f"s.add(If(x_{i}_{r}>0, y_{i}_{r}==1, y_{i}_{r}==0))")
    for i in range(c + 1, d1 + 1):
        for r in range(1, R + 1):
            _w(opfile, f"s.add(If(Ax_{i}_{r}>0, Ay_{i}_{r}==1, Ay_{i}_{r}==0))")
    for i in range(c + 1, d2 + 1):
        for r in range(1, R + 1):
            _w(opfile, f"s.add(If(Bx_{i}_{r}>0, By_{i}_{r}==1, By_{i}_{r}==0))")
    _w(opfile, "")


# =============================================================================
# MODULE 5 -- MIXER CONSISTENCY (ratio propagation, GLOBAL exponent 4^(i-1))
# =============================================================================

def mixerConsistency(c, d1, d2, R, opfile):
    """
    Uniform rule across stem AND both chains:
        Value_i_r = 4^(i-1) * (own fresh)_i_r + (incoming)_r

    Stem:
        R_1_r = x_1_r
        R_i_r = 4^(i-1) * x_i_r + t_i_(i-1)_r         i=2..c

    T1 (continues the SAME global exponent, references R_c_r directly):
        A_(c+1)_r = 4^c * Ax_(c+1)_r + tT_r           (tT_r = w1*R_c_r)
        A_i_r     = 4^(i-1) * Ax_i_r + At_i_(i-1)_r    i=c+2..d1

    T2 (same pattern, references the SAME R_c_r):
        B_(c+1)_r = 4^c * Bx_(c+1)_r + tU_r           (tU_r = w2*R_c_r)
        B_i_r     = 4^(i-1) * Bx_i_r + Bt_i_(i-1)_r    i=c+2..d2
    """
    _w(opfile, "# --- Stem ratio propagation ---")
    for r in range(1, R + 1):
        _w(opfile, f"s.add(R_1_{r}==x_1_{r})")
    for i in range(2, c + 1):
        for r in range(1, R + 1):
            _w(opfile,
               f"s.add(R_{i}_{r}==({N}**{i-1})*x_{i}_{r}+t_{i}_{i-1}_{r})")

    _w(opfile, "# --- T1 chain ratio propagation (references R_c directly) ---")
    for r in range(1, R + 1):
        _w(opfile,
           f"s.add(A_{c+1}_{r}==({N}**{c})*Ax_{c+1}_{r}+tT_{r})")
    for i in range(c + 2, d1 + 1):
        for r in range(1, R + 1):
            _w(opfile,
               f"s.add(A_{i}_{r}==({N}**{i-1})*Ax_{i}_{r}+At_{i}_{i-1}_{r})")

    _w(opfile, "# --- T2 chain ratio propagation (references SAME R_c) ---")
    for r in range(1, R + 1):
        _w(opfile,
           f"s.add(B_{c+1}_{r}==({N}**{c})*Bx_{c+1}_{r}+tU_{r})")
    for i in range(c + 2, d2 + 1):
        for r in range(1, R + 1):
            _w(opfile,
               f"s.add(B_{i}_{r}==({N}**{i-1})*Bx_{i}_{r}+Bt_{i}_{i-1}_{r})")
    _w(opfile, "")


# =============================================================================
# MODULE 6 -- TARGET VALIDITY (both in the same model)
# =============================================================================

def targetValidity(target1, target2, err, R, d1, d2, opfile):
    """
    S1 = sum(A_d1_j),  S2 = sum(B_d2_j)  -- actual denominators
    Cross-multiplied error tolerance (linear in Z3).
    """
    _w(opfile, "# --- T1 target validity ---")
    S1_expr = "+".join(f"A_{d1}_{r}" for r in range(1, R + 1))
    _w(opfile, f"S1 = {S1_expr}")
    _w(opfile, "s.add(S1 > 0)")
    for r in range(1, R + 1):
        _w(opfile, f"s.add({target1[r-1]}*S1-A_{d1}_{r}<={err}*S1)")
        _w(opfile, f"s.add(A_{d1}_{r}-{target1[r-1]}*S1<={err}*S1)")

    _w(opfile, "# --- T2 target validity ---")
    S2_expr = "+".join(f"B_{d2}_{r}" for r in range(1, R + 1))
    _w(opfile, f"S2 = {S2_expr}")
    _w(opfile, "s.add(S2 > 0)")
    for r in range(1, R + 1):
        _w(opfile, f"s.add({target2[r-1]}*S2-B_{d2}_{r}<={err}*S2)")
        _w(opfile, f"s.add(B_{d2}_{r}-{target2[r-1]}*S2<={err}*S2)")
    _w(opfile, "")


# =============================================================================
# MODULE 7 -- MINIMIZE
# =============================================================================

def minimizeObjective(c, d1, d2, R, opfile):
    """Minimize total active (reagent, mixer) pairs across stem + T1 + T2."""
    terms  = [f"y_{i}_{r}"  for i in range(1, c+1)     for r in range(1, R+1)]
    terms += [f"Ay_{i}_{r}" for i in range(c+1, d1+1)  for r in range(1, R+1)]
    terms += [f"By_{i}_{r}" for i in range(c+1, d2+1)  for r in range(1, R+1)]
    _w(opfile, f"s.minimize({'+'.join(terms)})")
    _w(opfile, "")


# =============================================================================
# ASSEMBLE: build the ONE combined SMT script
# =============================================================================

def buildCombinedSMT(c, d1, d2, target1, target2, err, R, z3opt, z3op):
    with open(z3opt, "w", encoding="utf-8") as opfile:
        initOPT(opfile)
        declareVariables(c, d1, d2, R, opfile)
        linearityConstraints(c, d1, d2, R, opfile)
        domainAndBudget(c, d1, d2, R, opfile)
        mixerConsistency(c, d1, d2, R, opfile)
        targetValidity(target1, target2, err, R, d1, d2, opfile)
        minimizeObjective(c, d1, d2, R, opfile)
        finishOPT(opfile, z3op)


# =============================================================================
# PRINT RESULT
# =============================================================================

def printResult(z3op, c, d1, d2, R):
    vals = _parseZ3(z3op)

    print(f"  Stem (leaf i=1 -> root/split i={c}):")
    for i in range(1, c + 1):
        xs = [vals.get(f"x_{i}_{r}", 0) for r in range(1, R+1)]
        Rs = [vals.get(f"R_{i}_{r}", 0) for r in range(1, R+1)]
        W  = vals.get(f"W_{i}_{i-1}", "-") if i > 1 else "leaf"
        print(f"    M{i}: x={xs}  W={W}  R={Rs}")

    w1 = vals.get('w1', '?'); w2 = vals.get('w2', '?')
    print(f"  Split: w1={w1} -> T1   w2={w2} -> T2")

    print(f"  T1 chain (leaf i={c+1} -> root i={d1}):")
    for i in range(c + 1, d1 + 1):
        xs = [vals.get(f"Ax_{i}_{r}", 0) for r in range(1, R+1)]
        As = [vals.get(f"A_{i}_{r}",  0) for r in range(1, R+1)]
        Wv = vals.get(f"AW_{i}_{i-1}", w1) if i > c+1 else w1
        print(f"    A{i}: Ax={xs}  W={Wv}  A={As}")
    S1 = sum(vals.get(f"A_{d1}_{r}", 0) for r in range(1, R+1))
    ratio1 = [round(vals.get(f"A_{d1}_{r}",0)/S1,4) if S1>0 else 0 for r in range(1,R+1)]
    print(f"  T1 ratio = {ratio1}  (S1={S1})")

    print(f"  T2 chain (leaf i={c+1} -> root i={d2}):")
    for i in range(c + 1, d2 + 1):
        xs = [vals.get(f"Bx_{i}_{r}", 0) for r in range(1, R+1)]
        Bs = [vals.get(f"B_{i}_{r}",  0) for r in range(1, R+1)]
        Wv = vals.get(f"BW_{i}_{i-1}", w2) if i > c+1 else w2
        print(f"    B{i}: Bx={xs}  W={Wv}  B={Bs}")
    S2 = sum(vals.get(f"B_{d2}_{r}", 0) for r in range(1, R+1))
    ratio2 = [round(vals.get(f"B_{d2}_{r}",0)/S2,4) if S2>0 else 0 for r in range(1,R+1)]
    print(f"  T2 ratio = {ratio2}  (S2={S2})")


# =============================================================================
# MAIN GENERATOR (combined single-SMT approach)
# =============================================================================

def generateSharedStemTreeSingleSMT(target1, target2, err,
                                     outDir="./sharedStem/",
                                     label="pair",
                                     maxDepth=13):
    """
    Find a shared-stem tree using ONE combined SMT per (c, d1, d2) attempt.
    R_c is a Z3 VARIABLE referenced by both T1 and T2 -- sharing is
    guaranteed by construction, no enumeration needed.
    """
    R = len(target1)
    assert len(target2) == R

    if all(abs(target1[j] - target2[j]) < 1e-9 for j in range(R)):
        print("T1 == T2: single tree satisfies both.")
        return {"note": "identical targets"}

    os.makedirs(outDir + "Z3Files", exist_ok=True)

    d1_ind = _findDepth(target1, err, R)
    d2_ind = _findDepth(target2, err, R)
    print(f"Individual depths: d1={d1_ind}  d2={d2_ind}")

    c_start = max(1, min(d1_ind, d2_ind) - 1)
    print(f"Starting c={c_start}, reducing to 1 on failure\n")

    for c in range(c_start, 0, -1):
        print(f"--- c={c} ---")
        for d1 in range(c + 1, maxDepth + 1):
            for d2 in range(c + 1, maxDepth + 1):
                print(f"  d1={d1} d2={d2}...", end=" ", flush=True)

                z3op  = os.path.join(outDir, "Z3Files",
                         f"{label}_c{c}_d1{d1}_d2{d2}_{err}").replace("\\", "/")
                z3opt = os.path.join(SCRIPT_DIR, "z3Combined.py").replace("\\", "/")

                buildCombinedSMT(c, d1, d2, target1, target2, err, R, z3opt, z3op)
                subprocess.call([sys.executable, z3opt])

                if _isSAT(z3op):
                    print("SAT OK\n")
                    printResult(z3op, c, d1, d2, R)
                    return {"c": c, "d1": d1, "d2": d2, "z3file": z3op}
                print("unsat")

        print()

    print("No shared stem found.")
    return None


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    result = generateSharedStemTreeSingleSMT(
        target1  = [0.25, 0.50, 0.25],
        target2  = [0.50, 0.25, 0.25],
        err      = 0.01,
        outDir   = "./sharedStem/",
        label    = "combined",
        maxDepth = 13,
    )
    if result and "note" not in result:
        print(f"\nDone: c={result['c']} d1={result['d1']} d2={result['d2']}")
