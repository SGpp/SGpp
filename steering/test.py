from pysgpp import *
import re
import tools

class Runner:
 
  dim = 3
  level = 5
  grid = None

  def __init__(self,fname):
    return

  ## build parable test function over [0,1]^d
  # @param dim dimension of the parable's space
  # @return returns a string that contains the function as string
  def buildParable(self,dim, timestep):
    function = "(x3+4)*(x3+4)*(x3+4)"
    
    for i in xrange(dim-1):
        function = function + "*x" + str(i+1) + "*(1-" + "x" + str(i+1) + ")"
    
    print "Generated Function: " + function    
    return function 


  ## evalutes a given function
  # @param function a string the gives the function; x1...xn must be the names of the placeholders
  # @param points sorted list of the coordinates (x1...xn) of evaluation point
  # @return returns the function value at points
  def evalFunction(self, function, points):
    for i in xrange(len(points)):
        function = re.sub("x" + str(i+1), points[i], function)
            
    return eval(function)

  def fillStorage(self):
    function = self.buildParable(self.dim, 4.0)
    self.grid = Grid.createLinearBoundaryGrid(self.dim)
    generator = self.grid.createGridGenerator()
    generator.regular(self.level)
    
    storage = self.grid.getStorage()
    node_values = DataVector(storage.size(), 1)

    
    for n in xrange(storage.size()):
      points = storage.get(n).getCoordsString().split()
      node_values[n] = self.evalFunction(function, points)
   
    alpha = DataVector(node_values);
    hierarchisation = self.grid.createOperationHierarchisation()
    hierarchisation.doHierarchisation(alpha)
    #hierarchisation.doDehierarchisation(alpha)    
    
    self.writeGnuplot("out.dat", self.grid, alpha, 100 )
    return
    
  def run(self):
    self.fillStorage()
    return

  def writeGnuplot(self, filename, grid, alpha, resolution):
    p = DataVector(1,3)
    fout = file(filename, "w")

    for z in xrange(100):
      for x in xrange(resolution):
        for y in xrange(resolution):
            p[0] = float(x) / (resolution - 1)
            p[1] = float(y) / (resolution - 1)
            p[2] = 0.01 * z
            pc = grid.createOperationEval().eval(alpha, p)
            fout.write("%f %f %f\n" % (p[0], p[1], pc))
    fout.close()
    return
# end Runner class

### Set up and start of the program
r = Runner("c")
r.run()
