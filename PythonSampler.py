#!/usr/bin/python

import subprocess
import shlex
import os
import time
import sys
import re

class AbstractParameter:
  PARAMETER_TYPES = ["define", "environment", "option"]
      
class CoresParameter(AbstractParameter):
  def __init__(self):
    self.min = 3
    self.max = 4
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
    if self.max == self.min:
      return "0.5"
    return str(float(self.value - self.min) / (self.max - self.min));
  
  def next(self):
    self.value += self.step

class BlocksizeParameter(AbstractParameter):
  def __init__(self):
    self.min = 64
    self.max = 256
    self.stepMult = 2
    self.value = self.min
    
  def __repr__(self):
    return "Blocksize"
  
  def getType(self):
    return "define"
  
  def reset(self):
    self.value = self.min
  
  def getVariableName(self):
    return "STREAMING_OCL_LOCAL_SIZE"
  
  def hasNext(self):
    # print self.value, "+", self.step, "->", (self.value + self.step)
    if self.value * self.stepMult > self.max:
      return False
    return True
  
  def getValue(self):
    return "STREAMING_OCL_LOCAL_SIZE=" + str(self.value)

  def getIndex(self):
    if self.max == self.min:
      return "0.5"
    return str(float(self.value - self.min) / (self.max - self.min));
  
  def next(self):
    self.value *= self.stepMult
    
class PrecisionParameter(AbstractParameter):
  
  def __init__(self):
    self.values = ["float", "double"]
    self.index = 0
    self.maxIndex = 1
  
  def __repr__(self):
    return "Precision"
  
  def getType(self):
    return "define"
  
  def reset(self):
    self.index = 0
    
  def hasNext(self):
    if self.index < self.maxIndex:
      return True
    return False
  
  def getIndex(self):
    if len(self.values) == 1:
      return "0.5"
    return str(float(self.index) / len(self.values));
  
  def getValue(self):
    return "STREAMING_OCL_INTERNAL_PRECISION=" + str(self.values[self.index])
  
  def next(self):
    self.index += 1

class UnrollMaxParameter(AbstractParameter):
  def __init__(self):
    self.min = 2
    self.max = 8
    self.step = 2
    self.value = self.min
    
  def __repr__(self):
    return "UnrollMax"
  
  def getType(self):
    return "define"
  
  def reset(self):
    self.value = self.min
  
  def getVariableName(self):
    return "STREAMING_OCL_MAX_DIM_UNROLL"
  
  def hasNext(self):
    if self.value + self.step > self.max:
      return False
    return True
  
  def getValue(self):
    return "STREAMING_OCL_MAX_DIM_UNROLL=" + str(self.value)

  def getIndex(self):
    if self.max == self.min:
      return "0.5"
    return str(float(self.value - self.min) / (self.max - self.min));
  
  def next(self):
    self.value += self.step

parameters = [
  CoresParameter(),
  BlocksizeParameter(),
  PrecisionParameter(),
  UnrollMaxParameter()
]

class Sampler:
  def __init__(self, resultFileName, parameters, jobs):
    self.resultFileName = resultFileName
    self.parameters = parameters
    self.jobs = str(jobs)

    self.bestParameters = None
    self.bestParameterDuration = sys.maxint
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
    my_env['LD_LIBRARY_PATH'] = 'lib/sgpp:lib/alglib'
    my_env['PYTHONPATH'] = 'lib/pysgpp'

    for envParameter in envParameters:
      keyValue = envParameter.getValue()
      my_env[keyValue[0]] = keyValue[1]
    
    # compile
    print "define parameter values:", " ".join([defineParameter.getValue() for defineParameter in defineParameters])
    command = "scons CPPFLAGS=\"" + ''.join([option.getValue() + " " for option in optionParameters]) + \
      "\" VERBOSE=1 -j" + self.jobs + " NO_UNIT_TESTS=1 CPPDEFINES=\"" + \
      " ".join([defineParameter.getValue() for defineParameter in defineParameters]) + "\" " + \
      "SG_PYTHON=0 SG_JAVA=0 OPT=1 CXX=g++-4.8 CC=g++-4.8"
    print command
    # command = ["./execTest.py"]
    print "origC:", command
    splittedCmd = shlex.split(command)    
    p = subprocess.Popen(splittedCmd, env=my_env) #stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
    p.wait()
    stdout, stderr = p.communicate()
    
    # execute
    command = ["bin/examples/datadriven/sampler"]
    startTimestamp = time.time()
    p = subprocess.Popen(command, shell=True, env=my_env)
    p.wait()
    stdout, stderr = p.communicate()
    duration = time.time() - startTimestamp
    print "duration:", duration

    if self.bestParameterDuration == None or duration < self.bestParameterDuration:
      self.bestParameterDuration = duration
      self.bestParameters = str([str(parameter) + "=" + str(parameter.getValue()) for parameter in optionParameters + envParameters + defineParameters])

    # extract runtime information
    dataTuple = [parameter.getIndex() for parameter in optionParameters + envParameters + defineParameters]
    print "configuration: ", str([str(parameter) + "=" + str(parameter.getValue()) for parameter in optionParameters + envParameters + defineParameters])
    #print "dataTuple:", dataTuple
    self.resultFile.write(", ".join(dataTuple) + ", " + str(duration) + "\n")
    
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

jobs = 1
for arg in sys.argv:
  match = re.search("JOBS=([0-9]+)", arg)
  if match != None:
    jobs = match.group(1)

print "number of jobs: ", jobs

sampler = Sampler(resultFileName, parameters, jobs)
sampler.sample()
