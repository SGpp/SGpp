from pysgpp import *
import re
import tools
import pickle

class Storage:
 
  grid = None
  node_value_vx = None
  node_value_vy = None
  node_value_vz = None
  
  alpha_vx = None
  alpha_vy = None  
  alpha_vz = None
  
  socketClient = None

  def __init__(self):
    self.deserializeGrid()
    self.hierarchise()

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

    fin = open("grid.vz", "r")
    value_vz = pickle.load(fin)

    storage = self.grid.getStorage()
    self.node_values_vx = DataVector(storage.size(), 1)
    self.node_values_vy = DataVector(storage.size(), 1)
    self.node_values_vz = DataVector(storage.size(), 1)

    for n in xrange(storage.size()):
      self.node_values_vx[n] = value_vx[n]
      self.node_values_vy[n] = value_vy[n]
      self.node_values_vz[n] = value_vz[n]

 
  def hierarchise(self):
    self.alpha_vx = DataVector(self.node_values_vx);
    self.alpha_vy = DataVector(self.node_values_vy);
    self.alpha_vz = DataVector(self.node_values_vz);

    hierarchisation = self.grid.createOperationHierarchisation()

    hierarchisation.doHierarchisation(self.alpha_vx)
    hierarchisation.doHierarchisation(self.alpha_vy)
    hierarchisation.doHierarchisation(self.alpha_vz)

    #hierarchisation.doDehierarchisation(alpha)    
    #self.writeGnuplot("out.dat", self.grid, self.alpha_vx, self.alpha_vy, self.alpha_vz, 100, [0.5,0.0])
    return

  def extractGrid(self, resolution, parameters):
    p = DataVector(1,5)
    values = []

    evaluator = self.grid.createOperationEval()
    
    for z in xrange(resolution):
      for x in xrange(resolution):
        for y in xrange(resolution):
          p[0] = float(x) / (resolution - 1)
          p[1] = float(y) / (resolution - 1)
          p[2] = float(z) / (resolution - 1)
          p[3] = parameters[0]
          p[4] = parameters[1]
          vx = evaluator.eval(self.alpha_vx, p)
          vy = evaluator.eval(self.alpha_vy, p)
          vz = evaluator.eval(self.alpha_vz, p)
          values.append(vx)
          values.append(vy)
          values.append(vz)

    return values
 
 
  def writeGnuplot(self, filename, grid, alpha_vx, alpha_vy, resolution, parameters):
    p = DataVector(1,4)
    fout = file(filename, "w")

    for y in xrange(resolution):
      for x in xrange(resolution):
        p[0] = float(x) / (resolution)
        p[1] = float(y) / (resolution)
        p[2] = parameters[0]
        p[3] = parameters[1]
        vx = grid.createOperationEval().eval(alpha_vx, p)
        vy = grid.createOperationEval().eval(alpha_vy, p)
        fout.write("%f %f %f %f\n" % (p[0], p[1], vx, vy))
    fout.close()
    return

# end Storage class

s = Storage()
