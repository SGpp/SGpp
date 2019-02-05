from combiScheme import *


def test():
    print("this is an original test by nico")
    return


def getstandardcombi(dim, b):
    print("Getting the Standard Combination level vectors and coefficients:")
    combi = CombiScheme(dim) 
    oldCombi = combi.getCombiScheme(1, b, dim)
    newCombi = []
    for x in oldCombi:
        a = [x[1]]
        a.extend(x[0])
        newCombi.append(a)
    return newCombi

