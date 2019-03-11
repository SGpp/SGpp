'''
Created on Aug 18, 2014

@author: franzefn
'''
from __future__ import division
from builtins import zip
from past.utils import old_div


def findEquivalent(sample, stats):
    # get sample in unit space
    p = tuple(sample.getExpandedUnit())
    found = p in stats
    if not found:
        min_diff = 10.
        # DAMN IT => this is shit!!!!
        keys = list(stats.keys())
        print( "search for equivalent for %s" % (p,) )
        j = 0
        while not found and j < len(keys):
            g = keys[j]
            diff1 = [old_div(abs(gi - pi), pi) for gi, pi in zip(g, p) if abs(pi) > 0]
            diff2 = [old_div(abs(gi - pi), gi) for gi, pi in zip(g, p) if abs(gi) > 0]
            diff = sum(diff1) + sum(diff2)
            min_diff = min(min_diff, diff)
            if diff < 1e-6:
                found = True
        if found:
            print( "found equivalent %s" % (g,) )
        else:
            print( "no equivalent found" )


import json
