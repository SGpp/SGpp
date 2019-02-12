from combiScheme import *


def test():
    print("this is an original test by nico")
    return

def newtestpart1():
    print("this is part 1 of another original test by nico")
    return CombiScheme(3)

def newtestpart2(schema):
    print("this is part 2 of another original test by nico")
    oldCombi = schema.getCombiScheme(1, 3, 3)
    newCombi = []
    for x in oldCombi:
        a = [x[1]]
        a.extend(x[0])
        newCombi.append(a)
    return newCombi

def getstandardcombi(dim, b):
    print("Space-SGDEAdapter: getstandardcombi")
    combi = CombiScheme(dim) 
    oldCombi = combi.getCombiScheme(1, b, dim)
    newCombi = []
    for x in oldCombi:
        a = [x[1]]
        a.extend(x[0])
        newCombi.append(a)
    return newCombi

def refineBlock(activeIndexSet):
    print("Space-SGDEAdapter: refineBlock")
