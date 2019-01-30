from combiScheme import *


def test():
    print("this is an original test by nico")
    return


def getstandardcombi(dim, b):
    combi = CombiScheme(dim)
    oldCombi = combi.getCombiScheme(1, b, dim)
    #print("That is the old data format")
    #print(oldCombi)
    newCombi = []
    for x in oldCombi:
        a = [x[1]]
        a.extend(x[0])
        newCombi.append(a)
    #print("That is the new data format")
    #print(newCombi)
    return newCombi

