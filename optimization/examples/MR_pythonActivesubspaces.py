import numpy as np
import matplotlib.pyplot as plt
import pysgpp

class ExampleFunction(pysgpp.OptScalarFunction):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(ExampleFunction, self).__init__(2)

    def eval(self, v):
        """Evaluates the function."""
        return np.exp(0.7*v[0] + 0.3*v[1])

pysgpp.OptPrinter.getInstance().setVerbosity(-1)
f = ExampleFunction()
gridType=pysgpp.GridType_NakBsplineModified
degree = 3
maxNumGridPointsMatrix = 100
initialLevel = 1

ASM = pysgpp.ASMatrixNakBspline(f,gridType,degree)
ASM.buildAdaptiveInterpolant(maxNumGridPointsMatrix,initialLevel,3)
ASM.createMatrixGauss()
ASM.evDecompositionForSymmetricMatrices()

eigenvalues = ASM.getEigenvaluesDataVector()
print(eigenvalues)
eigenvectors = ASM.getEigenvectorsDataMatrix()
print(eigenvectors)