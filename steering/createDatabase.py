from pysgpp import *
import re
import tools
import os, glob
import pickle

class Runner:
 
  dim = 3
  level = 7
  grid = None
  data = [ [ [[0 for bla in range(1)] for col in range(2**7) ] for row in range(2**7)] for bla2 in range(2**7)]
  node_value_vx = None
  node_value_vy = None

  def __init__(self,fname):
    return

  def readData(self):
    path = '../../DrivenCavity_128x128x128'
    steps = 0

    for infile in glob.glob( os.path.join(path, 'merged*data') ):
      for line in file(infile):
        line = line.split()
        x    = float(line[0])
        y    = float(line[1])
        vx   = float(line[3])
        vy   = float(line[4])
        
        self.data[int(x*(2**self.level-1))][int(y*(2**self.level-1))][steps] = [vx,vy]
      steps+=1
      print "processing " + infile

  def fillStorage(self):
    self.grid = Grid.createLinearBoundaryGrid(self.dim)
    generator = self.grid.createGridGenerator()
    generator.regular(self.level)
    
    storage = self.grid.getStorage()
    self.node_value_vx = [] #DataVector(storage.size(), 1)
    self.node_value_vy = [] #DataVector(storage.size(), 1)
    
    print "Creating regular grid with " + str(storage.size()) + " points."
 
    for n in xrange(storage.size()):
      points = storage.get(n).getCoordsString().split()
      self.node_value_vx.insert(n, self.data[int(float(points[0])*(2**self.level-1))][int(float(points[1])*(2**self.level-1))][int(float(points[2])*(2**self.level-1))][0])
      self.node_value_vy.insert(n, self.data[int(float(points[0])*(2**self.level-1))][int(float(points[1])*(2**self.level-1))][int(float(points[2])*(2**self.level-1))][1])
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
    return

  def hierarchise(self):
    alpha_vx = DataVector(self.node_value_vx2);
    alpha_vy = DataVector(self.node_value_vy2);
    hierarchisation = self.grid.createOperationHierarchisation()

    hierarchisation.doHierarchisation(alpha_vx)
    hierarchisation.doHierarchisation(alpha_vy)

    #hierarchisation.doDehierarchisation(alpha)    

    self.writeGnuplot("out2.dat", self.grid, alpha_vx, alpha_vy, 128)
    return

  def writeGnuplot(self, filename, grid, alpha_vx, alpha_vy, resolution):
    p = DataVector(1,3)
    fout = file(filename, "w")

    for z in xrange(100):
      for x in xrange(resolution):
        for y in xrange(resolution):
            p[0] = float(x) / (resolution - 1)
            p[1] = float(y) / (resolution - 1)
            p[2] = 0.01 * z
            vx = grid.createOperationEval().eval(alpha_vx, p)
            vy = grid.createOperationEval().eval(alpha_vy, p)
            fout.write("%f %f %f %f\n" % (p[0], p[1], vx, vy))
    fout.close()
    return

# end Runner class

### Set up and start of the program
r = Runner("c")
r.run()
