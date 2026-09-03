import graphviz
import LAFCADFL as ld

def getPlacement(W, d, pos):
    """
        Parameters:
        @Param1: Shared intermediate fluids between two node
        @Param2: Depth of the tree
        @Param3: Starting position of mixers
    """
    assignment = dict() # in the form of {M1: [timestamp, cells]}

    # At first place all the mixers on top of another (as skewed)
    for i in range(d-1, -1, -1):
        assignment[f"M{i}"] = [d-i-1, [[pos[0] , pos[1]  ],
                                      [pos[0]  , pos[1]+1],
                                      [pos[0]+1, pos[1]  ],
                                      [pos[0]+1, pos[1]+1]]
                              ]
    
    # Check wheather a non adjacent sharing exists or not.
    # If exists then move M{i-1} to desired cells based on units of reagents shared
    # i-1 as we are considering mixers starting from 1 to d
    for i in W:
        for j in W[i]:
            if i> j+1:
                if W[i][j] != 3:
                    for k in range(j+1, i): # Shift all the mixers in beteeen to the right
                        for p in assignment[f"M{k}"][1]:
                            p[1] += 1
                else:
                    for k in range(j+1, i): # Shift all the mixers in beteeen to the right-bottom
                        for p in assignment[f"M{k}"][1]:
                            p[0] += 1
                            p[1] += 1

    return assignment


def saveTree(x, W, ID, R=None, N=4):
    """
        Parameters:
        @Param1: Reagent usage per node (x)
        @Param2: Weight of each edges (W)
        @Param3: ID of the tree
        @Param4: Ratios per node (R), optional -- if given, each mixer node
                 is labeled with its actual computed ratio, not just its name
        @Param5: Mixer size N (default 4), used to normalize R into a
                 true fraction: node i's ratio = R[i][r] / N**(depth-i)
    """
    depth = len(x)
    # Create a graphviz object
    dot = graphviz.Digraph(comment=ID)
    dot.attr(rankdir='BT')

    # Add nodes for each mixer, labeled with its actual ratio if R is given
    for i in range(depth):
        if R is not None:
            denom = N ** (depth - i)
            ratios = [r / denom for r in R[i]]
            ratioStr = ":".join(f"{v:.3f}".rstrip('0').rstrip('.') or "0" for v in ratios)
            dot.node(f"M{i}", f"M{i}\n{ratioStr}")
        else:
            dot.node(f"M{i}")

    # Add edges from Mi to Mj
    for i in W:
        for j in W[i]:
            dot.edge(f"M{i}", f"M{j}", label=f"{W[i][j]}")

    # Add reagent units needed for each mixer
    for i, r in enumerate(x):
        for j, val in enumerate(r):
            if val != 0:
                dot.node(f"R{i+1}{j+1}", f"R{j+1}")
                dot.edge(f"R{i+1}{j+1}", f"M{i}", label=f"{val}")

    # return the object
    return dot


def isUnsatOrEmpty(z3File):
    """
    Returns True if the z3 output file does not contain a solved model
    (i.e. it is 'unsat', empty, or otherwise unparsable), False if it looks
    like a real 'name = value' model.
    """
    with open(z3File, 'r') as f:
        content = f.read().strip()
    if content == "" or content == "unsat":
        return True
    # A real model file's first non-blank line should contain " = "
    firstLine = content.splitlines()[0] if content else ""
    return " = " not in firstLine


def getMixerData(z3File, k, d):
    """
        Parameters:
        @Param1: z3 file name
        @Param2: Number of reagennts required
        @Param3: Depth of mixing tree
    """
    if isUnsatOrEmpty(z3File):
        raise ValueError(f"'{z3File}' does not contain a solved (SAT) model -- "
                          f"it is unsat, empty, or malformed. Check before calling getMixerData.")

    # R to store ratio at each intermediate node
    R = []
    # M to store what reagents and intermediate fluids are being used at each intermediate node
    M = dict()
    # x to store at each node units of reagent that are used
    x = []
    # W to store intermediate fluids that are shared between nodes
    W = dict()
    for i in range(d):
        R.append([0]*k)
        M[f"M{i}"] = []
        x.append([0]*k)

    # Open the z3file and store the values
    with open(z3File, 'r') as ipFile:
        for raw_line in ipFile:
            line = raw_line.rstrip('\n').split(" = ")

            # Skip blank lines, 'unsat', or any malformed line that isn't "name = value"
            if len(line) < 2:
                continue

            # Skip non-linearity helper variables (t_i_j_r)
            if line[0][0] == 't':
                continue

            val = int(line[1])
            RWx = line[0].split('_')

            # Skip anything that isn't in the expected "PREFIX_i_j" shape
            # (e.g. the standalone 'totalReagents' variable)
            if len(RWx) < 3:
                continue

            i, j = int(RWx[1])-1, int(RWx[2])-1

            if RWx[0] == 'R':
                R[i][j] = val
            elif RWx[0] == 'x':
                x[i][j] = val
                if val != 0:
                    M[f"M{i}"] += [f"R{j+1}"]*val
            elif RWx[0] == 'W':
                if val != 0:
                    if i not in W:
                        W[i] = dict()
                    W[i][j] = val
                    M[f"M{j}"] += [f"M{i}"]*val

    # Return the values
    return R, M, x, W


def saveTree_getArea(z3fileName, outputTreeImage, N, depth, startingCell):
    # R: mixer ratios
    # M: mixer components [R1, R2, R5, M1]
    # x: reagent usage
    # W: intermediate fluid sharing information
    R, M, x, W = getMixerData(z3fileName, N, depth)
    
    dot = saveTree(x, W, outputTreeImage, R=R, N=4)
    dot.render(outputTreeImage, format='png', cleanup=True)
    mixerData = getPlacement(W, depth, startingCell)

    assignment = dict()
    timestamp = dict()

    for key in mixerData:
        assignment[key] = mixerData[key][1]
        timestamp[mixerData[key][0]] = [key]

    # k, b, l = ld.KBL(assignment, M, timestamp)
    bb, area = 0, 0
    bb, area = ld.boundingbox(assignment)

    waste, reagentUsage = 0, 0
    # Calculate wastage
    reUse = dict()
    for i in W:
        if i not in reUse:
            reUse[i] = 0
        for j in W[i]:
            reUse[i] += W[i][j]

    for i in reUse:
        waste += 4-reUse[i]
    # Calculate total reagent usage
    for usage in x:
        reagentUsage += sum(usage)

    return waste, reagentUsage, bb, area


if __name__ == "__main__":
    id = 263
    z3path = f"./z3OutputFilesForShared/z3outputFile{id}"

    if isUnsatOrEmpty(z3path):
        print(f"id {id}: no SAT solution found in {z3path}, skipping")
    else:
        waste, reagentUsage, BB, area = saveTree_getArea(z3path, f"{id}", 4, 4, [4, 4])
        print(BB, area)
