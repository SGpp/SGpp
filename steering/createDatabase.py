from pysgpp import *
import re
import tools
import os, glob
import pickle
import numpy
from decimal import *


class Runner:
 
  dim = 5
  levels = [7,7,7,7,1]
  grid = None
  #data = [] # [ [ [ [[0 for bla in range(1) ] for col in range(2**levels[0])] for row in range(2**levels[1])] for bla2 in range(2**levels[2])] for bla2 in range(2**levels[3])]
  data_x = numpy.zeros([2**levels[4]+1,2**levels[1]+1,2**levels[2]+1,2**levels[3]+1,2**levels[0]+1] ,float)
  data_y = numpy.zeros([2**levels[4]+1,2**levels[1]+1,2**levels[2]+1,2**levels[3]+1,2**levels[0]+1] ,float)
  node_value_vx = None
  node_value_vy = None

  def __init__(self,fname):
    return

  def readData(self):
    path = '/work/butnaru/SparseData/output'
    steps = 0
    Re = 0

    print "Will create data store of size ", 2**self.levels[0]+1, 2**self.levels[1]+1, 2**self.levels[2]+1, 2**self.levels[3]+1

    data = glob.glob( os.path.join(path, 'merged*data') )
    data.sort()

    for infile in data:
      print "processing " + infile, "step " , steps, " Re ", Re
      if steps == 129:
        steps = 0
        Re+=1
 	
      for line in file(infile):
        line = line.split()
        x    = float(line[0])
        y    = float(line[1])
        z    = float(line[2])

        vx   = float(line[3])
        vy   = float(line[4])
      
        ix = int(round(x*(2**self.levels[0])))
        iy = int(round(y*(2**self.levels[1])))
        iz = int(round(z*(2**self.levels[2])))
        #print "index", ix, iy, x, y, steps, Re, vx, vy
  
        self.data_x[Re][steps][ix][iy][iz]  = vx
        self.data_y[Re][steps][ix][iy][iz]  = vy

      steps+=1
  
  def fillStorage(self):
    self.grid = Grid.createLinearTrapezoidBoundaryGrid(self.dim)
    generator = self.grid.createGridGenerator()
    generator.regularWithDifferentLevels(self.levels)
    
    storage = self.grid.getStorage()
    self.node_value_vx = []
    self.node_value_vy = []
    
    print "Creating regular grid with " + str(storage.size()) + " points."
 
    for n in xrange(storage.size()):
      points = storage.get(n).getCoordsString().split()
      i0 = int(float(points[0])*(2**self.levels[0]))
      i1 = int(float(points[1])*(2**self.levels[1]))
      i2 = int(float(points[2])*(2**self.levels[2]))
      i3 = int(float(points[3])*(2**self.levels[3]))
      i4 = int(float(points[4])*(2**self.levels[4]))
 
      self.node_value_vx.insert(n, self.data_x[i4][i3][i0][i1][i2])
      self.node_value_vy.insert(n, self.data_y[i4][i3][i0][i1][i2])
    return
 
  def serializeGrid(self):
    text = self.grid.serialize()
    fout = open("grid.obstacle.store", "w")
    fout.write(text)
    fout.close()
    
    fout = open("grid.obstacle.vx", "w")
    pickle.dump(self.node_value_vx, fout)

    fout = open("grid.obstacle.vy", "w")
    pickle.dump(self.node_value_vy, fout)

  def run(self):
    print "Reading data."
    self.readData()
    print "Filling storage."
    self.fillStorage()
    print "Serializing storage."
    self.serializeGrid()

    print "Finished"
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
