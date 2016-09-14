'''
Created on May 18, 2016

@author: franzefn
'''
import numpy as np
from matplotlib2tikz import save as tikz_save

from pysgpp import Grid, SLinearBase
from pysgpp.extensions.datadriven.uq import mpl, plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize, \
    isHierarchicalAncestor, isHierarchicalAncestorDimx
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d, plotGrid2d


def plotBasis(basis, l, i, n=1000, w=1, **kws):
    x = np.linspace(0, 1, n)
    y = [w * basis.eval(l, i, xi) for xi in x]
    plt.plot(x, y, **kws)


def evalInequality(grid, basis, x, gppp, gpkk, gpll=None):
    gpp, pl, pi, px = gppp
    gpk, kl, ki, kx = gpkk
    if gpll is None:
        # 1d case
        return basis.eval(pl, pi, x) - basis.eval(pl, pi, kx) * basis.eval(kl, ki, x)
    else:
        # 2d case
        gpl, ll, li, lx = gpll

        res = basis.eval(pl, pi, x) - basis.eval(pl, pi, kx) * basis.eval(kl, ki, x)

        if isHierarchicalAncestor(gpk, gpl):
            res -= (basis.eval(pl, pi, lx) - basis.eval(pl, pi, kx) * basis.eval(kl, ki, lx)) * basis.eval(ll, li, x)
        else:
            res -= basis.eval(pl, pi, lx) * basis.eval(ll, li, x)

        return res

# plot results for numDims = 1
numDims = 1
level = 4

grid = Grid.createLinearGrid(numDims)
grid.getGenerator().regular(level)
gs = grid.getStorage()

basis = SLinearBase()

p = 0
gpp = gs.getPoint(p)
idim = 0
pl, pi, px = gpp.getLevel(idim), gpp.getIndex(idim), gpp.getStandardCoordinate(idim)

mpl.rcParams['axes.labelsize'] = 35

for weighted in [True, False]:
    fig = plt.figure()
    ax = plt.subplot(111)
    plotBasis(basis, pl, pi, label=r"$\phi_%i(x)$" % p)

    for l in xrange(3, 7):
        gpp, gpl = gs.getPoint(p), gs.getPoint(l)
        ll, li, lx = gpl.getLevel(idim), gpl.getIndex(idim), gpl.getStandardCoordinate(idim)

        if weighted:
            w = basis.eval(pl, pi, lx)
            label = r"$\phi_%i(x_%i) \phi_%i(x)$" % (p, l, l)
        else:
            w = 1
            label = r"$\phi_%i(x)$" % (l)

        plotBasis(basis, ll, li, w=w, label=label)

        x = np.linspace(0, 1, 1000)
        y = [evalInequality(grid, basis, xi, (gpp, pl, pi, px), (gpl, ll, li, lx))
             for xi in x]
        # plt.plot(x, y, label="res")

        # plt.vlines([px, lx], 0, 1)

    # Place a legend to the right of this smaller subplot.
    plt.legend()

    if weighted:
        tikz_save("src/positive_weighted_%i_%id.tex" % (idim, numDims),
                  encoding="utf8")
    else:
        tikz_save("src/positive_unweighted_%i_%id.tex" % (idim, numDims),
                  encoding="utf8")
    fig.show()

# # plot results for numDims = 2
# numDims = 2
# level = 4
#
# grid = Grid.createLinearGrid(numDims)
# grid.getGenerator().regular(level)
# gs = grid.getStorage()
#
# p, k, l = 1, 9, 43
#
# fig = plt.figure()
# plotGrid2d(grid, show_numbers=True)
# plt.title("p, k, l = (%i, %i, %i)" % (p, k, l))
# plt.savefig("figures/positive_reference_grid_%id.pdf" % numDims)
# fig.show()
#
# gpp, gpk, gpl = gs.getPoint(p), gs.getPoint(k), gs.getPoint(l)
#
# for idim in xrange(numDims):
#     pl, pi, px = gpp.getLevel(idim), gpp.getIndex(idim), gpp.getStandardCoordinate(idim)
#     kl, ki, kx = gpk.getLevel(idim), gpk.getIndex(idim), gpk.getStandardCoordinate(idim)
#     ll, li, lx = gpl.getLevel(idim), gpl.getIndex(idim), gpl.getStandardCoordinate(idim)
#
#     basis = SLinearBase()
#
#     fig = plt.figure()
#     plotBasis(pl, pi, label=r"$\phi_i(x)$")
#     plotBasis(kl, ki, w=basis.eval(pl, pi, kx), label=r"$\phi_i(x_k) \phi_k(x)$")
#     plotBasis(ll, li, w=basis.eval(pl, pi, lx), label=r"$\phi_i(x_l) \phi_l(x)$")
#
#     x = np.linspace(0, 1, 1000)
#     y = [evalInequality(grid, basis, xi, (gpp, pl, pi, px), (gpk, kl, ki, kx), (gpl, ll, li, lx))
#          for xi in x]
#     plt.plot(x, y, label="res")
#
#     plt.vlines([px, kx, lx], 0, 1)
#
#     plt.title("d = %i: is ancestor? %s" % (idim, isHierarchicalAncestorDimx(gpk.getLevel(idim),
#                                                                             gpk.getIndex(idim),
#                                                                             gpl.getLevel(idim),
#                                                                             gpl.getIndex(idim))))
#     plt.legend()
#     plt.savefig("figures/positive_%i_%id.pdf" % (idim, numDims))
#     fig.show()
#
plt.show()
