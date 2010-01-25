from pysgpp import *
import re
import tools
import pickle
from client2 import SocketClient

class Runner:
 
  dim = 2
  level = 7
  grid = None
  node_value_vx = None
  node_value_vy = None
  
  alpha_vx = None
  alpha_vy = None  
  
  socketClient = None

  def __init__(self,fname):
    return

  def hierarchise(self):
    self.alpha_vx = DataVector(self.node_values_vx);
    self.alpha_vy = DataVector(self.node_values_vy);
    hierarchisation = self.grid.createOperationHierarchisation()

    hierarchisation.doHierarchisation(self.alpha_vx)
    hierarchisation.doHierarchisation(self.alpha_vy)

    #hierarchisation.doDehierarchisation(alpha)    
    #self.writeGnuplot("out.dat", self.grid, alpha_vx, alpha_vy, 128)
    return

  def deserializeGrid(self):
    fin = open("grid.store", "r")
    text = fin.read()
    self.grid = Grid.unserialize(text)
    fin.close()

    fin = open("grid.vx", "r")
    value_vx = pickle.load(fin)

    fin = open("grid.vy", "r")
    value_vy = pickle.load(fin)

    storage = self.grid.getStorage()
    self.node_values_vx = DataVector(storage.size(), 1)
    self.node_values_vy = DataVector(storage.size(), 1)

    for n in xrange(storage.size()):
      self.node_values_vx[n] = value_vx[n]
      self.node_values_vy[n] = value_vy[n]
 
  def run(self):
    self.deserializeGrid()
    self.hierarchise()
    self.socketClient = SocketClient()
    counter = 0 

    while 1:
      #params = self.socketClient.recv_end()
      #print params
      data = self.extractGrid(100, counter * 0.01)
      self.socketClient.sendCurrentGrid(data)

      counter = counter + 1

      if counter > 100:
        counter = 0

    return

  def extractGrid(self, resolution, parameterValue):
    p = DataVector(1,3)
    gridAsString = ""
    
    for x in xrange(resolution):
      for y in xrange(resolution):
        p[0] = float(x) / (resolution - 1)
        p[1] = float(y) / (resolution - 1)
        p[2] = parameterValue
        vx = self.grid.createOperationEval().eval(self.alpha_vx, p)
        vy = self.grid.createOperationEval().eval(self.alpha_vy, p)
        gridAsString = gridAsString + " " +  str(vx)
        gridAsString = gridAsString + " " +  str(vy)
    return gridAsString
  
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
