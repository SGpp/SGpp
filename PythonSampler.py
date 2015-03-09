import subprocess
import shlex
import os
import time

class AbstractParameter:
  PARAMETER_TYPES = ["define", "environment", "option"]
    
class VectorParameter(AbstractParameter):
  def __init__(self):
    self.values = ["-msse3", "-mavx", ""]
    self.index = 0
  
  def __repr__(self):
    return "Vector"
  
  def getType(self):
    return "option"
  
  def reset(self):
    self.index = 0
    
  def hasNext(self):
    if self.index < 2:
      return True
    return False
  
  def getIndex(self):
    return str(float(self.index) / len(self.values));
  
  def getValue(self):
    return self.values[self.index]
  
  def next(self):
    self.index += 1
      
class CoresParameter(AbstractParameter):
  def __init__(self):
    self.min = 1
    self.max = 3
    self.step = 1
    self.value = self.min
  
  def __repr__(self):
    return "Cores"
  
  def getType(self):
    return "environment"
  
  def reset(self):
    self.value = self.min
  
  def hasNext(self):
    if self.value + self.step > self.max:
      return False
    return True
  
  def getValue(self):
    return ("OMP_NUM_THREADS", str(self.value))
  
  def getIndex(self):
    return str(float(self.value - self.min) / (self.max - self.min));
  
  def next(self):
    self.value += self.step

class BlocksizeParameter(AbstractParameter):
  def __init__(self):
    self.min = 1
    self.max = 3
    self.step = 1
    self.value = self.min
    
  def __repr__(self):
    return "Blocksize"
  
  def getType(self):
    return "define"
  
  def reset(self):
    self.value = self.min
  
  def getVariableName(self):
    return "BLOCKSIZE"
  
  def hasNext(self):
    # print self.value, "+", self.step, "->", (self.value + self.step)
    if self.value + self.step > self.max:
      return False
    return True
  
  def getValue(self):
    return "BLOCKSIZE=" + str(self.value)

  def getIndex(self):
    return str(float(self.value - self.min) / (self.max - self.min));
  
  def next(self):
    self.value += 1

parameters = [
  VectorParameter(),
  CoresParameter(),
  BlocksizeParameter()
]

class Sampler:
  def __init__(self, resultFileName, parameters):
    self.resultFileName = resultFileName
    self.parameters = parameters
  def dimIter(self):
    index = 0
    while index < len(parameters):
      parameter = parameters[index]
      if parameter.hasNext():
        values = [paraIter.getValue() for paraIter in self.parameters]
        yield parameters
        
        parameter.next()
        
        for i in range(0, index):
          parameters[i].reset()
        index = 0
      else:
        index += 1
    yield parameters
    
  def execute(self, optionParameters, envParameters, defineParameters):
    # print optionParameters
    # print envParameters
    # print defineParameters
    
    # setup environment
    my_env = os.environ.copy()
    for envParameter in envParameters:
      keyValue = envParameter.getValue()
      my_env[keyValue[0]] = keyValue[1]
    
    # compile
    command = "scons CPPFLAGS=\"" + ''.join([option.getValue() + " " for option in optionParameters]) + \
      "\" VERBOSE=1 -j4 NO_UNIT_TESTS=1 CPPDEFINES=\"" + \
      "".join([defineParameter.getValue() for defineParameter in defineParameters]) + "\""
    # command = ["./execTest.py"]
    print "origC:", command
    splittedCmd = shlex.split(command)
    print "split:", splittedCmd
    p = subprocess.Popen(splittedCmd, env=my_env) #stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
    stdout, stderr = p.communicate()
    
    # execute
    command = ["datadriven/examples/sampler"]
    startTimestamp = time.time()
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=my_env)
    stdout, stderr = p.communicate()
    duration = time.time() - startTimestamp
    print "duration:", duration

    # extract runtime information
    dataTuple = [parameter.getIndex() for parameter in optionParameters + envParameters + defineParameters]
    print "dataTuple:", dataTuple
    self.resultFile.write(", ".join(dataTuple) + ", " + str(duration))
    
  def sample(self):
    self.resultFile = open(self.resultFileName, "w")
    for parameters in self.dimIter():
      # print "parameters:", parameters
      
      optionParameters = []
      envParameters = []
      defineParameters = []
      
      for parameter in parameters:
        if parameter.getType() == "option":
          optionParameters += [parameter]
        elif parameter.getType() == "environment":
          envParameters += [parameter]
        elif parameter.getType() == "define":
          defineParameters += [parameter]
        else:
          print "error"
          
      self.execute(optionParameters, envParameters, defineParameters)
      
    self.resultFile.close()

resultFileName = "samples.dat"

sampler = Sampler(resultFileName, parameters)
sampler.sample()
