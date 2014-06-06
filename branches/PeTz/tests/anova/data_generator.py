#!/usr/bin/python
# -*- coding: utf-8 -*-
''' Generator of artificial data for ANOVA decomposition test '''


from numpy import cos, sin, pi, exp

if __name__ == '__main__':
    f1 = open('refinement_strategy_sum.csv', 'wt')
    f2 = open('refinement_strategy_x1.csv', 'wt')
    f3 = open('refinement_strategy_sincos.csv', 'wt')

    print >> f1, 'x,y,target'
    print >> f2, 'x,y,target'
    print >> f3, 'x,y,target'

    for x in [float(x) / 100 for x in xrange(101)]:
        for y in [float(y) / 100 for y in xrange(101)]:
            f1.write('%.3f, %.3f, %.8f\n' % (x, y, x ** 2 + y ** 2))
            f2.write('%.3f, %.3f, %.8f\n' % (x, y, exp(2 * x ** 2) + y))
            f3.write('%.3f, %.3f, %.8f\n' % (x, y, sin((x + .5) * 2 * pi)
                     * cos((y + .5) * 2 * pi)))
    f1.close()
    f2.close()
    f3.close()

