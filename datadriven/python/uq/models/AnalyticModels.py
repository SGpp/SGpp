from Model import Model
import numpy as np
from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder

class Parabola(Model):
    
    def buildUQSetting(self):
        f = lambda x, **kws: np.prod([4 * xi * (1 - xi) for xi in x])
        self.uqBuilder.withSimulation(f)

class Atan(Model):

    def buildUQSetting(self):
        f = lambda x, **kws: np.atan(50 * (x[0] - .35)) + pi / 2 + 4 * x[1] ** 3 + exp(x[0] * x[1] - 1)
        self.uqBuilder.withSimulation(f)
