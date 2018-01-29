import numpy as np


def findSetBits(x):
    # find which bits are set in an interger
    ans = []
    bx = bin(x)[2:]
    n = len(bx)
    for i, b in enumerate(bx[::-1]):
        if b == '1':
            ans.append(i)
    return tuple(ans)


def sortPermutations(perms, index_return=False):
    """
    Sort perms with respect (1) to their length and (2) their lexical order
    @param perms:
    """
    ans = [None] * len(perms)
    indices = np.ndarray(len(perms), dtype="int")
    ix = 0
    for n in np.sort(np.unique([len(key) for key in perms])):
        # load subset of perms with length n
        nperms = {}
        for i, perm in enumerate(perms):
            if len(perm) == n:
                tperm = tuple(perm)
                nperms[tperm] = i

        for perm in sorted(nperms.keys()):
            ans[ix] = perm
            indices[ix] = nperms[perm]
            ix += 1

    if index_return:
        return ans, indices
    else:
        return ans

def computeTotalEffects(sobol_indices):
    total_effects = {}
    for k in [perm for perm in sobol_indices.keys() if len(perm) == 1]:
        total_effects[k] = 0.0
        for perm, sobol_index in sobol_indices.items():
            if k[0] in perm:
                total_effects[k] += sobol_index
    return total_effects
