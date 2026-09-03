try:
    from .LAFCA import *
    from .DFL import *
    from .NTM import *
except ImportError:
    # Fallback for when this file is run directly (no parent package context)
    # instead of imported through LAFCADFL as a package.
    from LAFCA import *
    from DFL import *
    from NTM import *

import itertools

maxDeg = -1
blockagePlacement = []
intermediateFluid = []



def dfs(cell, cells, placement, units, blockage, index, blockedPos, fluidPos, grid=None):
    if units[index] == 0:
        if index+1 == len(cells):
            global maxDeg, blockagePlacement, intermediateFluid
            internalDeg = 0
            externalDeg = 0
            d = [-1, 0, 1, 0, -1]
            row = len(grid) if grid is not None else 0
            col = len(grid[0]) if grid is not None and row > 0 else 0
            for [x, y] in blockedPos:
                for i in range(4):
                    nx, ny = x+d[i], y+d[i+1]
                    if [nx, ny] in blockedPos:
                        internalDeg += 1
                    elif (grid is not None) and (0 <= nx < row) and (0 <= ny < col) \
                            and (grid[nx][ny] != '*'):
                        # Reward touching an already-placed external blockage
                        # (e.g. a non-adjacent-shared M cell from an earlier
                        # timestep). Without this, every equal-size subset of
                        # a 2x2 box scores identically on internal adjacency
                        # alone, so the free/leftover cell ends up chosen by
                        # iteration order rather than by which side is safe --
                        # this is what left the fresh reagent's cell adjacent
                        # to a persisted blockage wall with no clean route.
                        externalDeg += 1
            deg = (internalDeg // 2) + externalDeg
            if maxDeg< deg:
                maxDeg=deg
                blockagePlacement = blockedPos[:]
                intermediateFluid = fluidPos[:]
            return
        dfs(cells[index+1], cells, placement, units, blockage, index+1, blockedPos, fluidPos, grid)
    
    dir = [0, 0, -1, 0, 1, 0]
    for i in range(5):
        x, y = cell[0]+dir[i], cell[1]+dir[i+1]
        if ([x, y] in placement[index]) & ([x, y] not in blockedPos):
            blockedPos.append([x, y])
            units[index]-=1
            fluidPos.append(blockage[index])
            dfs([x, y], cells, placement, units, blockage, index, blockedPos, fluidPos, grid)
            units[index]+=1
            fluidPos.pop()
            blockedPos.pop()



def findBlockages(loadingCells, units, blockages, grid=None):
    for element in itertools.product(*loadingCells):
        dfs(element[0], element, loadingCells, units, blockages, 0, [], [], grid)



def getPlacementAndLoading(Mixtures, parentMix, loadingCells, reagentList, blockageList, units, grid):
    '''
        Intermediate fluids are considered as blockages
    '''
    # Make all the blockings adjacent
    '''
        Structure of blockageList is:
        blockageList["M1"] = {"M4": {"M1": [[1,2], [1,3], [2,2]], "M2": [[2,2]]} },
                    {"M5": {"M3": [[2,4]]} },
                    {"M6": {}
        }
        units["M4"] = [
            {"M1":2, "M2":1},
            {"M3":1},
            {},
        ]
    '''
    # If the intermediate fluid required and allocated space is equal then there's no choice but put
    # the intermediate fluids in that allocated space only.
    for _ in range(0, 2):
        for mixture in loadingCells:
            removeMix = []
            for mix in blockageList[mixture]:
                if units[mixture][mix] == len(blockageList[mixture][mix]):
                    removeMix.append(mix)
                    for x, y in blockageList[mixture][mix]:
                        # Fix the intermediate cell pos in grid and remove it from loading cells
                        grid[x][y] = mix
                        loadingCells[mixture].remove([x, y])
                        # remove [x, y] in other blockages also
                        for mix1 in blockageList[mixture]:
                            if (mix != mix1) & ([x, y] in blockageList[mixture][mix1]):
                                blockageList[mixture][mix1].remove([x, y])
            for mix in removeMix:
                del blockageList[mixture][mix]
                del units[mixture][mix]

    # Generate all possible combination and choose best out of it
    ind = 0
    parent = dict()
    blockageLoad = []
    blockageUnits = []
    blockageNames = []
    for mix in blockageList:
        ind+=1
        for child in blockageList[mix]:
            parent[child] = mix 
            blockageLoad.append(blockageList[mix][child][:])
            blockageUnits.append(units[mix][child])
            blockageNames.append(child)

    global maxDeg, blockagePlacement, intermediateFluid
    if len(blockageNames)> 0:
        findBlockages(blockageLoad, blockageUnits, blockageNames, grid)
        # print(blockagePlacement)
        # print(intermediateFluid)

    for i in range(len(blockagePlacement)):
        loadingCells[parent[intermediateFluid[i]]].remove(blockagePlacement[i])
        grid[blockagePlacement[i][0]][blockagePlacement[i][1]] = intermediateFluid[i]

    maxDeg = -1
    blockagePlacement = []
    intermediateFluid = []

    # All the cells that need to load
    cellsToLoad = []
    for mixture in loadingCells:
        for cell in loadingCells[mixture]:
            cellsToLoad.append(cell)

    # List of all the reagents
    allReagents = []
    for mix in reagentList:
        for reagent in reagentList[mix]:
            allReagents.append(reagent)

    # Check wheather z3 call is necessary or not , i.e if only one cell is empty then no need 
    # or all reagents are same
    toDel = []
    for mixture in reagentList:
        reagents = set()
        for reagent in reagentList[mixture]:
            reagents.add(reagent)
        if len(reagents) <= 1:
            if len(reagents) == 1:
                for cell in loadingCells[mixture]:
                    grid[cell[0]][cell[1]] = reagentList[mixture][0]
            toDel.append(mixture)
    
    for mixture in toDel:
        del reagentList[mixture]
        del loadingCells[mixture]

    toDel = []
    for mixture in loadingCells:
        if len(loadingCells[mixture]) == 0:
            toDel.append(mixture)
    
    for mix in toDel:
        del loadingCells[mixture]

    # Used LAFCA to place the reagents in each cell
    i = 0
    for mix in loadingCells:
        assignment = createFile(reagentList[mix], loadingCells[mix], 'z3File.py', 'output'+str(i)+'.txt')
        i += 1
        for reagent in assignment:
            for j in range(len(assignment[reagent])):
                # Updating the grid
                x = assignment[reagent][j][0]
                y = assignment[reagent][j][1]
                grid[x][y] = reagent

    for r in grid:
        print(r)
    print()
    print(Mixtures)
    # Need to make row and col as local variable that can be passed in DFL
    loadingPaths = DFL([0,9],[9,0],grid,allReagents,cellsToLoad)
    # print(loadingPaths)
    for mix in Mixtures:
        for x, y in Mixtures[mix]:
            if [x, y] in parentMix[mix]:
                grid[x][y] = mix
            else:
                grid[x][y] = '*'   # washed: no longer needed by any future consumer

    for r in grid:
        print(r)
    print()

    global v
    v = -1

    return loadingPaths



###################################
def getMix(root):
    '''
        @param1 NTM root
        
        @returns: a dictionary that contains all the mixtures and their ratio list
        Simple BFS implementation
    '''
    queue = deque()
    queue.append(root)
    mixture = dict()

    while queue:
        s = len(queue)
        nodes = []
        while s:
            s -= 1
            node = queue.popleft()
            reags = []
            nodes += [child for child in node.children if child.children != []]
            for child in node.children:
                reags.extend([child.value]*child.volume)
            mixture[node.value] = reags

        for n in nodes:
            queue.append(n)
    return mixture



def boundingbox(assignments):
    """
        To find the bounding box by finding the Xmax, Xmin, Ymax, Ymin
        @param1: Assignments of all the mixers from NTM output.
        @return: Pair of integers, {Total cells required, Area of the boundinng box}
    """
    uniqueCells = set()
    x_max, y_max, x_min, y_min = 0, 0, 15, 15
    for mix in assignments:
        for cell in assignments[mix]:
            uniqueCells.add(''.join(str(coord) for coord in cell))
            x_max = max(x_max, cell[0])
            x_min = min(x_min, cell[0])
            y_max = max(y_max, cell[1])
            y_min = min(y_min, cell[1])
    
    area = (x_max - x_min + 1) * (y_max - y_min + 1)
    return len(uniqueCells), area
            


def KBL(assignment, mixtures, timestamp):
    '''
        Returns KBL information of a mixing tree
    '''
    allFlow, allBendings, allLengths = 0, 0, 0
    row, col = 10, 10
    grid = []
    parent = dict() # Stores parent coordinates
    for mix in mixtures:
        for reagent in set(mixtures[mix]):
            if reagent[0] == 'M':
                i, j = int(reagent[1]), int(mix[1])
                cnt = 0
                for rgnt in mixtures[mix]:
                    cnt += (rgnt == reagent)

                if reagent not in parent:
                    if i> j+1:
                        if cnt < 3:
                            parent[reagent] = [assignment[mix][0], assignment[mix][2]]
                        else:
                            parent[reagent] = [assignment[mix][0], assignment[mix][1], assignment[mix][2]]
                    else:
                        parent[reagent] = assignment[mix][:]
                else:
                    if i> j+1:
                        if cnt < 3:
                            parent[reagent] += [assignment[mix][0], assignment[mix][2]]
                        else:
                            parent[reagent] += [assignment[mix][0], assignment[mix][1], assignment[mix][2]]
                    else:
                        parent[reagent] += assignment[mix][:]

    parent['M0'] = assignment['M0'][:]

    for _ in range(row):
        grid.append(['*']*col)
        
    # Get the parallel loading cells in each time stamp
    for t in timestamp:
        # Make the grid
        Mixtures = dict()
        loadingCells = dict() # Cells that participate in mixing at current timestamp
        reagentList = dict() # List of reagents used in mixture M1
        blockageList = dict() # List of intermediate fluids in each mixture
        units = dict() # Units of intermediate fluids needed in mixtures

        print("timestamp", t)
        for mix in timestamp[t]: # Mi Mj etc
            Mixtures[mix] = assignment[mix][:]
            loadingCells[mix] = assignment[mix][:]
            reagents = []
            blockage = dict() # Store blockage list for each mixture Mi
            unit = dict() # Store units of intermediate fluids required
            for reagent in mixtures[mix]:
                if reagent[0] == 'M': # indicates intermediate fluid (blockage)
                    if reagent not in blockage:
                        i, j = int(reagent[1]), int(mix[1])
                        unit[reagent] = 0
                        for rgnt in mixtures[mix]:
                            unit[reagent] += (rgnt == reagent)
                        
                        # [top-left, top-right, bottom-left, bottom-right]
                        blockage[reagent] = []
                        
                        if i> j+1:
                            if unit[reagent] == 1:
                                blockage[reagent].append(assignment[reagent][2])
                            elif unit[reagent] == 2:
                                blockage[reagent].append(assignment[reagent][0])
                                blockage[reagent].append(assignment[reagent][2])
                            elif unit[reagent] == 3:
                                blockage[reagent].append(assignment[reagent][0])
                                blockage[reagent].append(assignment[reagent][1])
                                blockage[reagent].append(assignment[reagent][2])
                        else:
                            for cell in assignment[reagent]:
                                if cell in assignment[mix]:
                                    blockage[reagent].append(cell)
                                # else:
                                #     grid[cell[0]][cell[1]] = '*' # Washing
                else:
                    reagents.append(reagent)
            # Reagents are in reagents
            reagentList[mix] = reagents #list
            # Blockages are in blockage and their positions
            blockageList[mix] = blockage #dict
            # Intermediate fluid units are stored in units
            units[mix] = unit #dict

            # if t == 4:
            #     print(reagentList[mix])
            #     print(blockageList[mix])
            #     print(units[mix])

        loadingPaths = getPlacementAndLoading(Mixtures, parent, loadingCells, reagentList, blockageList, units, grid)
        totalPathLength, totalBendings = 0, 0
        for order in loadingPaths:
            totalBendings += order[1]
            totalPathLength += len(order[2])
            print(order[0], 'Bendings:', order[1], 'Path Length:', len(order[2]))
        
        print('Flow:', len(loadingPaths), ',Total Bendings:', totalBendings, ',Total Path Length:', totalPathLength)
        allFlow += len(loadingPaths)
        allBendings += totalBendings
        allLengths += totalPathLength
    
    print('K', allFlow, 'B', allBendings, 'L', allLengths)
    return allFlow, allBendings, allLengths



def getPlacementAndTimestamp(root):
    '''
        Input: Root is a list format mentioned in NTM package from which NTM root is generated.
        Generate the tree from the list provided
        Use NTM to get the placement of the tree and time stamp at which each mixture will execute
    '''
    ntmroot = listToTree(root)
    output_assignment_set = ntm(ntmroot, [5, 5], [1]) # returns [moduleID, timeStamp, Binding, WashSequence] for every sequence

    # Get the corrospondence mixture reagents and intermediate fluids
    mixture = getMix(ntmroot)

    # Assignment of all the internal node
    assignment = {}
    # timestamp at which particular mixture is going to execute
    timeStamp = {}
    for item in output_assignment_set:
        if item[0][0] == 'M':
            assignment[item[0]] = item[2]
            if item[1] not in timeStamp:
                timeStamp[item[1]] = [item[0]] # i.e timeStamp[1] = ["M1", "M2", "M5"]
            else:
                timeStamp[item[1]].append(item[0])

    BB, area = boundingbox(assignment)
    print("area ", area)
    K, B, L = KBL(assignment, mixture, timeStamp)
    return BB, area, K, B, L
