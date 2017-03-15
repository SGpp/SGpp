from Model import Model

class Parabola(Model):
    
    def evaluate(self, sample):
        return np.prod([4 * xi * (1 - xi) for xi in sample.getActiveUnit()])

