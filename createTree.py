from NTM import *
row=10
col=10
def createNode(rootList, x, W, ind):
    """
        Create a rootList from the list data of reagents in each level and sharing of intermediate
        fluid data.
    """
    if ind == len(x):
        return [], 0, 0, 0
    rootList.append(W[ind])
    for i, k in enumerate(x[ind]):
        if k != 0:
            rootList+=[f'R{i+1}']*k
    children = []
    createNode(children, x, W, ind+1)
    if len(children) != 0:
        rootList.append(children)
    return rootList

def getKBL(x, W, ind):
    # As the tree is skewed so only specific numbers of permutation of fluids are possible in a mixer
    cache = dict()
    cache['1 1 1 1'] = [4, 9, 4*(row+col-1)]
    cache['1 1 2'] = [3, 6, 3*(row+col-1)]
    cache['1 3'] = [2, 5, 2*(row+col-1)]
    cache['2 2'] = [2, 4, 2*(row+col-1)]
    cache['4'] = [1, 5, row+col+1]
    cache['1 1 1'] = [3, 7, 3*(row+col-1)]
    cache['1 2'] = [2, 4, 2*(row+col-1)]
    cache['3'] = [1, 3, 1*(row+col-1)]
    cache['1 1'] = [2, 5, 2*(row+col-1)]
    cache['2'] = [1, 2, 1*(row+col-1)]
    cache['1'] = [1, 3, 1*(row+col-1)]

    if ind == len(x):
        return 0, 0, 0
    K, B, L = 0, 0, 0
    child = []
    for i, k in enumerate(x[ind]):
        if k != 0:
            child.append(k)
    child = sorted(child)
    print(child)
    key = ' '.join(str(val) for val in child)
    K += cache[key][0]
    B += cache[key][1]
    L += cache[key][2]
    k, b, l = getKBL(x, W, ind+1)
    K += k
    B += b
    L += l
    return K, B, L


def createRoot(fileName, depth, r, f=0):
    R = []
    x = []
    W = [4]*(depth)
    for _ in range(depth):
        R.append([0]*r)
        x.append([0]*r)

    with open(fileName, "r+") as fp:
        lines = fp.read().split('\n')
        for line in lines[:-1]:
            var, val = line.split('=')
            vars = var.split('_')
            if vars[0] == 'R':
                R[int(vars[1])-1][int(vars[2])-1] = int(val)
            elif vars[0] == 'x':
                x[int(vars[1])-1][int(vars[2])-1] = int(val)
            else:
                W[int(vars[1])-1] = int(val)

    """
        In case one want the KBL information of a skewed mixing tree f=1 and we use predetermined values to calculate
    """
    K, B, L = 0, 0, 0
    if f==1:
        K, B, L =  getKBL(x, W, 0)

    # Code to get list of the tree
    rootList = []
    createNode(rootList, x, W, 0)
    # Creating list to tree that is already implemented in NTM library
    print(rootList)
    root = listToTree(rootList)
    """
        Return K,B,L along with root if needed
    """
    return K, B, L, root if f==1 else root


def showTree(fileName, outputFileName, depth, r):
    root = createRoot(fileName, depth, r)
    saveTree(root, f'./OutputTree/{outputFileName}.png')

def getDepth(z3file):
    depth = 1
    mixingRatioLen = 0
    with open(z3file, 'r') as fp:
        line = fp.readline()
        while line:
            if line[0] == 'W':
                line = line.split(" = ")
                pref = line[0].split("_")
                i, j = int(pref[1]), int(pref[2])
                if i == j+1:
                    depth+=1
            elif line[0] == 'R':
                line = line.split(' = ')
                line = line[0].split('_')
                mixingRatioLen = max(mixingRatioLen, int(line[2]))
            line = fp.readline()
    
    return depth, mixingRatioLen


###########################################################################

def createSkewedTree(z3File, r, N, outputFileName):
    depth, mixingRatioLen = getDepth(z3File)




if __name__ == "__main__":
    getDepth('z3outputFile')