from combiScheme import *

def initadaptivescheme(dim, level):
    scheme = CombiScheme(dim)
    scheme.init_adaptive_combi_scheme(level, 1)
    return scheme

def getcombischeme(scheme):
    oldCombi =  scheme.getCombiScheme(1,1,1,do_print=False)
    return formatCombi(oldCombi)

def isrefinable(scheme, levelvec):
    #print("levelvec: ",levelvec)
    return scheme.is_refinable(sgdecombitolevelvec(levelvec))

def refineblock(scheme, levelvec):
    print("Space-SGDEAdapter: refineblock")
    scheme.update_adaptive_combi(sgdecombitolevelvec(levelvec))
    return scheme

def formatCombi(oldCombi):
    newCombi = []
    for x in oldCombi:
        a = [x[1]]
        a.extend(x[0])
        newCombi.append(a)
    return newCombi

def sgdecombitolevelvec(sgdecombi):
    sgdecombi.pop(0)
    return sgdecombi
