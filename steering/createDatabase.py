from pysgpp import *
import re
import tools
import os, glob
import pickle
import numpy

class Runner:
 
  dim = 4
  levels = [7,7,7,3]
  grid = None
  data = [] # [ [ [ [[0 for bla in range(1)] for col in range(2**levels[0])] for row in range(2**levels[1])] for bla2 in range(2**levels[2])] for bla2 in range(2**levels[3])]
  data_x = numpy.zeros([2**levels[0],2**levels[1],2**levels[2],2**levels[3]],float)
  data_y = numpy.zeros([2**levels[0],2**levels[1],2**levels[2],2**levels[3]],float)
  node_value_vx = None
  node_value_vy = None

  def __init__(self,fname):
    return

  def readData(self):
    path = '../../DrivenCavity_128x128x128'
    steps = 0
    Re = 0

    for infile in glob.glob( os.path.join(path, 'merged*data') ):
      if steps < 128:    	
        for line in file(infile):
          line = line.split()
          x    = float(line[0])
          y    = float(line[1])
          vx   = float(line[3])
          vy   = float(line[4])
      
          #print "index", int(x*(2**self.levels[0]-1)), int(y*(2**self.levels[1]-1)), steps, Re 
          self.data_x[int(x*(2**self.levels[0]-1))][int(y*(2**self.levels[1]-1))][steps][Re] = vx
          self.data_y[int(x*(2**self.levels[0]-1))][int(y*(2**self.levels[1]-1))][steps][Re] = vy
        steps+=1
        print "processing " + infile
      else:
        steps = 0
        Re+=1

  def fillStorage(self):
    self.grid = Grid.createLinearBoundaryGrid(self.dim)
    generator = self.grid.createGridGenerator()
    generator.regularWithDifferentLevels(self.levels)
    
    storage = self.grid.getStorage()
    self.node_value_vx = []
    self.node_value_vy = []
    
    print "Creating regular grid with " + str(storage.size()) + " points."
 
    for n in xrange(storage.size()):
      points = storage.get(n).getCoordsString().split()
      self.node_value_vx.insert(n, self.data_x[int(float(points[0])*(2**self.levels[0]-1))][int(float(points[1])*(2**self.levels[1]-1))][int(float(points[2])*(2**self.levels[2]-1))][int(float(points[3])*(2**self.levels[3]-1))])
      self.node_value_vy.insert(n, self.data_y[int(float(points[0])*(2**self.levels[0]-1))][int(float(points[1])*(2**self.levels[1]-1))][int(float(points[2])*(2**self.levels[2]-1))][int(float(points[3])*(2**self.levels[3]-1))])
    return
 
  def serializeGrid(self):
    text = self.grid.serialize()
    fout = open("grid.store", "w")
    fout.write(text)
    fout.close()
    
    fout = open("grid.vx", "w")
    pickle.dump(self.node_value_vx, fout)

    fout = open("grid.vy", "w")
    pickle.dump(self.node_value_vy, fout)

  def run(self):
    print "Reading data."
    self.readData()
    print "Filling storage."
    self.fillStorage()
    print "Serializing storage."
    self.serializeGrid()

    self.hierarchise()
    return

  def hierarchise(self):
    alpha_vx = DataVector(1, len(self.node_value_vx))
    alpha_vy = DataVector(1, len(self.node_value_vy))
    
    for i in range(len(self.node_value_vx)):
      alpha_vx[i] = self.node_value_vx[i] 
      alpha_vy[i] = self.node_value_vy[i] 
      
    hierarchisation = self.grid.createOperationHierarchisation()

    hierarchisation.doHierarchisation(alpha_vx)
    hierarchisation.doHierarchisation(alpha_vy)

    #hierarchisation.doDehierarchisation(alpha)    

    self.writeGnuplot("out.dat", self.grid, alpha_vx, alpha_vy, 128)
    return

  def writeGnuplot(self, filename, grid, alpha_vx, alpha_vy, resolution):
    p = DataVector(1,4)
    fout = file(filename, "w")

    for z in xrange(100):
      for x in xrange(resolution):
        for y in xrange(resolution):
            p[0] = float(x) / (resolution - 1)
            p[1] = float(y) / (resolution - 1)
            p[2] = 1
            p[3] = 0.01 * z
            vx = grid.createOperationEval().eval(alpha_vx, p)
            vy = grid.createOperationEval().eval(alpha_vy, p)
            fout.write("%f %f %f %f\n" % (p[0], p[1], vx, vy))
    fout.close()
    return

# end Runner class

### Set up and start of the program
r = Runner("c")
p = DataVector(1,4)
print p
r.run()
