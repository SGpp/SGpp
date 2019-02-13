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

def initadaptivescheme(dim, level):
    print("Space-SGDEAdapter: initadaptivscheme")
    scheme = CombiScheme(dim)
    scheme.init_adaptive_combi_scheme(level, 1)
    return scheme

def getcombischeme(scheme):
    print("Space-SGDEAdapter: getcombischeme")
    oldCombi =  scheme.getCombiScheme(1,1,1)
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

def refineBlock(currentScheme):
    
    print("Space-SGDEAdapter: refineBlock")
