#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import sys

import numpy as np
import pysgpp

try:
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D
  skipPlot = False
except ImportError:
  skipPlot = True

def getGridPoints(gridStorage):
  N, dim = gridStorage.getSize(), gridStorage.getDimension()
  L = np.zeros((N, dim), dtype=np.uint64)
  I = np.zeros((N, dim), dtype=np.uint64)

  for k in range(N):
    for d in range(dim):
      L[k,d] = gridStorage.getPoint(k).getLevel(d)
      I[k,d] = gridStorage.getPoint(k).getIndex(d)

  return I / 2**L

dim, n, p, hasBoundary = 2, 4, 3, True
f = lambda XX: (((20*XX[:,1]-5)-5.1*(20*XX[:,0]-5)**2/(4*np.pi**2)+5*(20*XX[:,0]-5)/np.pi-6)**2 +
                10*(1-1/(8*np.pi))*np.cos(20*XX[:,0]-5) + 10)

basis1d = pysgpp.SBsplineBase(p)
basis = pysgpp.HeterogeneousBasis(dim, basis1d)
combinationGrid = pysgpp.CombinationGrid.fromRegular(dim, n, basis, hasBoundary)
gridStorage = pysgpp.HashGridStorage(dim)
combinationGrid.combinePoints(gridStorage)
gridPointsSG = getGridPoints(gridStorage)
valuesSparseGrid = pysgpp.DataVector(f(gridPointsSG))
valuesFullGrids = pysgpp.DataVectorVector()
combinationGrid.distributeValuesToFullGrids(gridStorage, valuesSparseGrid, valuesFullGrids)

operationPole = pysgpp.OperationPoleVector()
pysgpp.OperationPoleHierarchisationGeneral.fromHeterogenerousBasis(basis, operationPole)

for doFullGrid in [True, False]:
  if doFullGrid:
    fullGridIndex = 1
    fullGrid = combinationGrid.getFullGrids()[fullGridIndex]
    level = fullGrid.getLevel()
    operationHierarchisation = pysgpp.OperationUPFullGrid(fullGrid, operationPole)
    valuesFullGrid = valuesFullGrids[fullGridIndex]
    surplusesFullGrid = pysgpp.DataVector(valuesFullGrid)
    operationHierarchisation.apply(surplusesFullGrid)
  else:
    operationHierarchisation = pysgpp.OperationUPCombinationGrid(combinationGrid, operationPole)
    surplusesFullGrids = pysgpp.DataVectorVector(valuesFullGrids)
    operationHierarchisation.apply(surplusesFullGrids)

  xx0 = np.linspace(0, 1, 65)
  xx1 = np.linspace(0, 1, 65)
  XX0, XX1 = np.meshgrid(xx0, xx1)
  XX = np.column_stack([XX0.flatten(), XX1.flatten()])
  XXdm = pysgpp.DataMatrix(XX)
  YYdv = pysgpp.DataVector(0)

  if doFullGrid:
    operationEvaluation = pysgpp.OperationEvalFullGrid(fullGrid)
    operationEvaluation.eval(surplusesFullGrid, XXdm, YYdv)
  else:
    operationEvaluation = pysgpp.OperationEvalCombinationGrid(combinationGrid)
    operationEvaluation.eval(surplusesFullGrids, XXdm, YYdv)

  YY = np.array([YYdv[k] for k in range(YYdv.getSize())])
  YY = YY.reshape(XX0.shape)

  if skipPlot:
    print("Skipping plot due to failed import of Matplotlib.")
    sys.exit(0)

  if doFullGrid:
    Xdm = pysgpp.DataMatrix(0, 0)
    pysgpp.IndexVectorRange(fullGrid).getPoints(fullGrid, Xdm)
    X = np.array([[Xdm.get(k, j) for j in range(Xdm.getNcols())] for k in range(Xdm.getNrows())])
  else:
    X = gridPointsSG

  fig = plt.figure(figsize=(6, 6))
  ax = fig.gca(projection="3d")
  ax.plot_surface(XX0, XX1, YY)
  ax.plot(X[:,0], X[:,1], "k.", zs=f(X), ms=10)

plt.show()
