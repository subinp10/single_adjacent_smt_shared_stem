import createTreeForShared as cts
import createTree as ct


def KBL(cache, mix):
    hashedKey = ' '.join(str(val) for val in sorted(mix) if val> 0)
    return cache[hashedKey]


def KBL_0(row, col, mix):
    # For leaf node
    cache = dict()
    cache['1 1 1 1'] = [4, 9, 4*(row+col-1)]
    cache['1 1 2'] = [3, 6, 3*(row+col-1)]
    cache['1 3'] = [2, 5, 2*(row+col-1)]
    cache['2 2'] = [2, 4, 2*(row+col-1)]
    cache['4'] = [1, 5, row+col+1]

    return KBL(cache, mix)

    
def KBL_1(row, col, mix):
    # For intermediate nodes that have 1 unit of shared intermediate fluids
    cache = dict()
    cache['1 1 1'] = [3, 7, 3*(row+col-1)]
    cache['1 2'] = [2, 4, 2*(row+col-1)]
    cache['3'] = [1, 3, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_2(row, col, mix):
    # For intermediate nodes that have 2 units of shared intermediate fluids
    cache = dict()
    cache['1 1'] = [2, 5, 2*(row+col-1)]
    cache['2'] = [1, 2, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_3(row, col, mix):
    # For intermediate nodes that have 3 units of shared intermediate fluids
    cache = dict()
    cache['1'] = [1, 3, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_0_1(row, col, mix):
    # For creation of nodes that need to create on side cells with 1 unit of shared fluid
    # Non-neighbouring sharing 2 units
    cache = dict()
    cache['1 1 1'] = [3, 8, 3*(row+col-1)]
    cache['1 2'] = [2, 5, 2*(row+col-1)]
    cache['3'] = [1, 4, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_0_2(row, col, mix):
    # For creation of nodes that need to create on side cells with 2 unit of shared fluid
    # Non-neighbouring sharing 2 units
    cache = dict()
    cache['1 1'] = [2, 5, 2*(row+col-1)]
    cache['2'] = [1, 2, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_0_3(row, col, mix):
    # For creation of nodes that need to create on side cells with 3 unit of shared fluid
    # Non-neighbouring sharing 1-2 units
    cache = dict()
    cache['1'] = [1, 3, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_1_1(row, col, mix):
    # For creation of nodes that need to create on side cells with 1 unit of shared fluid
    # Non-neighbouring sharing 3 units
    cache = dict()
    cache['1 1 1'] = [3, 7, 3*(row+col-1)]
    cache['1 2'] = [2, 4, 2*(row+col-1)]
    cache['3'] = [1, 3, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_1_2(row, col, mix):
    # For creation of nodes that need to create on side cells with 22 unit of shared fluid
    # Non-neighbouring sharing 3 units
    cache = dict()
    cache['1 1'] = [2, 5, 2*(row+col-1)]
    cache['2'] = [1, 2, 1*(row+col-1)]

    return KBL(cache, mix)


def KBL_1_3(row, col, mix):
    # For creation of nodes that need to create on side cells with 3 unit of shared fluid
    # Non-neighbouring sharing 3 units
    cache = dict()
    cache['1'] = [1, 3, 1*(row+col-1)]

    return KBL(cache, mix)


def getLoadingData(fileName, coord, d, row, col, N, k):
    if N != 4:
        print("Only works for mixer-4")
        return -1, -1, -1

    K, B, L = 0, 0, 0 # Return values
    cx, cy = coord

    R, M, x, W = cts.getMixerData(fileName, k, d)
    assignment = cts.getPlacement(W, d, [cx, cy])
    print(R)
    print(M)

    timeStamp = dict()
    for key in assignment:
        scheduleTime = assignment[key][0]
        if scheduleTime not in timeStamp:
            timeStamp[scheduleTime] = []
        timeStamp[scheduleTime].append(key)

    for t in range(d):
        for mixer in timeStamp[t]:
            coordinates = assignment[mixer][1]
            print(mixer)
            print(coordinates)

            intermediateFluids = 0

            for m in M[mixer]:
                if m[0] == 'M':
                    intermediateFluids += 1

            if intermediateFluids == N:
                continue

            callKBL = ""
            if (coordinates[0][0]==cx) & (coordinates[0][1]==cy):
                callKBL = f'KBL_{intermediateFluids}'
            elif (coordinates[0][0]==cx) & (coordinates[0][1]==cy+1):
                callKBL = f'KBL_0_{intermediateFluids}'
            elif (coordinates[0][0]==cx+1) & (coordinates[0][1]==cy+1):
                callKBL = f'KBL_1_{intermediateFluids}'
            if callKBL in globals():
                res = globals()[callKBL](row, col, x[int(mixer[1:])])
                K += res[0]
                B += res[1]
                L += res[2]
            else:
                print("Function not found")
    return K, B, L

if __name__ == "__main__":
    id = 1014
    depth = 5
    fileName = f"./z3OutputFilesForShared/z3outputFile{id}"
    print(getLoadingData(fileName, [5, 5], depth, 10, 10, 4, 5))
