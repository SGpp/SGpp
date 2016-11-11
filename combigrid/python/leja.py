# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:40:46 2015

@author: Julian
"""

import numpy as np

def maximum_leja(next_point, z, lower_bound, upper_bound, weight=lambda x : 1):
    """
    the maximums function for the leja points
    the next leja point is the input so that this function returns
    its maximum
    """
    if (next_point[0] < lower_bound or next_point[0] > upper_bound):
        return -100
    prod = 1
    for i in range(len(z)):
        prod *= np.abs(next_point[0] - z[i])
    return prod * weight(next_point[0])

def invert_maximum_leja(next_point, z, lower_bound, upper_bound, weight=lambda x : 1):
    """
    The Nelder-Mead Simplex algorithm used to find the next leja point searches
    for the minimum of the function, so we just invert our maximum function
    """
    return (-1) * maximum_leja(next_point, z, lower_bound, upper_bound, weight)

def leja_poly(next_point, z, lower_bound, upper_bound, weight=lambda x : 1):
    if (next_point[0] < lower_bound or next_point[0] > upper_bound):
        return 100
    prod = 1
    for i in range(len(z)):
        prod *= (next_point[0] - z[i])
    return weight(next_point[0]) * prod

def leja_points(start, count, lower_bound, upper_bound, weight=lambda x : 1, debug=False, view=False):
    """
    calculates the next COUNT leja points with START = z_0
    returns the leja points in a list
    """
    leja = [start]
    if view:
        print "Leja Points:"
    for p in range(count):
        # calculate next leja point
        tries = []
        values = []

        if debug:
            print "-----------------------"

        for i in range(len(leja)):
            leja = sorted(leja)
            low = lower_bound
            if (i > 0):
                low = leja[i-1]
            up = leja[i]

            if (debug):
                print "Untere Grenze: " + str(low)
                print "Obere Grenze: " + str(up)

            f = lambda x : invert_maximum_leja(x, leja, low, up, weight)

            tries.extend(calc_min(f, low, up))
            values.append(f([tries[-1]]))

            if debug:
                print "Try: " + str(tries[-1])
                print "Value: " + str(values[-1])

        if not leja[-1] == upper_bound:
            if (debug):
                print "Untere Grenze: " + str(leja[-1])
                print "Obere Grenze: " + str(upper_bound)
            f = lambda x : invert_maximum_leja(x, leja, leja[-1], upper_bound, weight)
            tries.extend(calc_min(f, leja[-1], upper_bound))
            values.append(f([tries[-1]]))
            if debug:
                print "Try: " + str(tries[-1])
                print "Value: " + str(values[-1])

        leja.append(tries[values.index(min(values))])

        if debug:
            print "Taken: " + str(leja[-1])
        if view:
            if not p == count - 1:
                print str(leja[-1]) + ", "
            else:
                print str(leja[-1])
    return sorted(leja)

def calc_min(f, lower_bound, upper_bound):
    startvalue = (upper_bound - lower_bound) / 2 + lower_bound
    return minimize(f, startvalue, method='nelder-mead').x

if __name__ == '__main__':
    weight = lambda x : np.sin(x * np.pi)
    start = 0.5
    lower_bound = 0
    upper_bound = 1
    count = 8
    try:
        from scipy.optimize import minimize
        points = leja_points(start, count, lower_bound, upper_bound, weight)
        print "Sorted:"
        print points
    except:
        pass
