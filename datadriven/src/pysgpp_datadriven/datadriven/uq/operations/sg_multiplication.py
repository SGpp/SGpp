#
#
# def mult(grid, alpha, f, epsilon=1e-10, refnums=0, pointsNum=10):
#     # copy grid
#     n_grid = copyGrid(grid, level=8, deg=1)
#     n_gs = n_grid.getStorage()
#     n_gridGen = n_grid.createGridGenerator()
#
#     basis_alpha = DataVector(alpha)
#
#     # dehierarchization
#     n_gs = n_grid.getStorage()
#     ps = DataMatrix(n_gs.size(), n_gs.dim())
#     p = DataVector(n_gs.dim())
#     for i in xrange(n_gs.size()):
#         n_gs.get(i).getCoords(p)
#         ps.setRow(i, p)
#
#     nodalValues = evalSGFunctionMulti(grid, basis_alpha, ps)
#     for i in xrange(n_gs.size()):
#         ps.getRow(i, p)
#         nodalValues[i] = f(p.array(), nodalValues[i])
#
#     n_alpha = hierarchize(n_grid, nodalValues)
#
#     # TODO: Improve this refinement
#     # now refine adaptively refnum times
#     for _ in xrange(refnums):
#         rp = max(1, min(pointsNum, n_gridGen.getNumberOfRefinablePoints()))
#         n_gridGen.refine(SurplusRefinementFunctor(n_alpha, rp, epsilon))
#         # extend alpha vector...
#         basis_alpha.resizeZero(n_gs.size())
#         # ... initialize it with function values and ...
#         nodalValues = evalGridPoints(n_grid, basis_alpha, f)
#         # ... do hierarchization
#         n_alpha = hierarchize(n_grid, nodalValues)
#
#         # # -----------------------------------------
#         # # check if interpolation property is given
#         # p = DataVector(n_gs.dim())
#         # for i in xrange(n_gs.size()):
#         #     gp = n_gs.get(i)
#         #     gp.getCoords(p)
#         #     y1 = nodalValues[i]
#         #     y2 = evalSGFunction(n_grid, n_alpha, p)
#         #     if abs(y1 - y2) >= 1e-13:
#         #         print "i:", p, y1, y2, abs(y1 - y2)
#         #         # assert abs(y1 - y2) < 1e-12
#         # # -----------------------------------------
#
#         if len(n_alpha) == n_gs.size():
#             break
#
#     return n_grid, n_alpha
