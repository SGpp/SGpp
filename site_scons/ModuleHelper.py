# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import fnmatch
import inspect
import os

class Module(object):
  def __init__(self, moduleDependencies, additionalDependencies=[],
               additionalBoostTestDependencies=[], excludeFiles=[]):
    self.moduleDependencies = moduleDependencies
    self.additionalDependencies = additionalDependencies
    self.additionalBoostTestDependencies = additionalBoostTestDependencies
    self.excludeFiles = excludeFiles
    self.cpps = []
    self.hpps = []
    self.objs = []

    variables = inspect.stack()[1][0].f_globals
    for name, value in variables.iteritems():
      globals()[name] = value

    env.Export("moduleDependencies")
  
  def scanSource(self, sourceFolder="src"):
    sourceFolder = os.path.join(Dir(".").abspath, sourceFolder)

    for currentFolder, subdirNames, fileNames in os.walk(sourceFolder, topdown=True):
      subdirNames.sort()
      fileNames.sort()
      sconscriptPath = os.path.join(currentFolder, "SConscript")
      if (currentFolder != sourceFolder) and os.path.exists(sconscriptPath):
        module = self
        env.SConscript(sconscriptPath, exports="module")
        # remove subfolders from iteration as they are already processed
        # (this is why topdown=True is required)
        subdirNames[:] = []
      else:
        for fileName in fnmatch.filter(fileNames, "*.cpp"):
          if fileName in self.excludeFiles: continue
          cpp = os.path.join(currentFolder, fileName)
          self.cpps.append(cpp)
          self.objs.append(env.SharedObject(cpp))
        for fileName in fnmatch.filter(fileNames, "*.hpp"):
          if fileName in self.excludeFiles: continue
          hpp = os.path.join(currentFolder, fileName)
          self.hpps.append(hpp)

    for hpp in self.hpps:
      headerSourceList.append(os.path.join(moduleName, hpp))
      headerDestList.append(hpp.split(os.sep, 1)[1])

  def buildLibrary(self):
    libname = "sgpp%s" % moduleName
    if env["USE_STATICLIB"]:
      libname += "static"
      self.moduleDependencies = [module + "static" for module in self.moduleDependencies]
    env.Export("libname")
    self.libname = libname

    if env["USE_STATICLIB"]:
      libsuffix = env["LIBSUFFIX"]
      self.lib = env.StaticLibrary(target=self.libname, source=self.objs, LIBPATH=BUILD_DIR,
                                   LIBS=self.moduleDependencies + self.additionalDependencies)
    else:
      libsuffix = env["SHLIBSUFFIX"]
      self.lib = env.SharedLibrary(target=self.libname, source=self.objs, LIBPATH=BUILD_DIR,
                                   LIBS=self.moduleDependencies + self.additionalDependencies)

    for module in self.moduleDependencies:
      if module.startswith("sgpp"):
        otherLib = os.path.join("#", BUILD_DIR.path, env["LIBPREFIX"] + module + libsuffix)
        env.Depends(self.lib, otherLib)

    self.libInstall = env.Install(BUILD_DIR, self.lib)
    libraryTargetList.append(self.lib)
    installTargetList.append(self.libInstall)

  def buildExamples(self, exampleFolder="examples", additionalExampleDependencies=[]):
    exampleEnv = env.Clone()
    exampleEnv.AppendUnique(LIBS=[self.libname] + self.moduleDependencies +
                                  self.additionalDependencies + additionalExampleDependencies)
    for fileName in os.listdir(exampleFolder):
      if fnmatch.fnmatch(fileName, "*.cpp"):
        cpp = os.path.join(exampleFolder, fileName)
        self.cpps.append(cpp)
        example = exampleEnv.Program(source=cpp)
        exampleEnv.Depends(example, self.libInstall)
        exampleTargetList.append(example)
      elif fnmatch.fnmatch(fileName, "*.hpp"):
        hpp = os.path.join(exampleFolder, fileName)
        self.hpps.append(hpp)

  def runPythonTests(self):
    if not env["NO_UNIT_TESTS"] and env["SG_PYTHON"]:
      moduleTest = env.Test(os.path.join("tests", "test_{}.py".format(moduleName)))
      testTargetList.append(moduleTest)

  def buildBoostTests(self, boostTestFolder="tests", compileFlag="COMPILE_BOOST_TESTS"):
    if env[compileFlag]:
      testEnv = env.Clone()
      testEnv.AppendUnique(LIBS=self.moduleDependencies + self.additionalDependencies +
                                 [self.libname] + ["boost_unit_test_framework"] +
                                 self.additionalBoostTestDependencies)

      testObjs = []
      for currentFolder, subdirNames, fileNames in os.walk(boostTestFolder, topdown=True):
        for fileName in fnmatch.filter(fileNames, "*.cpp"):
          cpp = os.path.join(currentFolder, fileName)
          self.cpps.append(cpp)
          testObjs.append(testEnv.SharedObject(cpp))
        for fileName in fnmatch.filter(fileNames, "*.hpp"):
          hpp = os.path.join(currentFolder, fileName)
          self.hpps.append(hpp)
      if len(testObjs) > 0:
        self.boostTestExecutable = \
            os.path.join(boostTestFolder, "test_{}_boost".format(moduleName)) + \
            (".exe" if env["PLATFORM"] == "win32" else "")
        test = testEnv.Program(self.boostTestExecutable, testObjs)
        testEnv.Depends(test, self.libInstall)

  def runBoostTests(self, boostTestFolder="tests",
                    compileFlag="COMPILE_BOOST_TESTS", runFlag="RUN_BOOST_TESTS"):
    if env[compileFlag] and env[runFlag]:
      testRun = env.BoostTest(self.boostTestExecutable + "_run", source=self.boostTestExecutable)
      boostTestTargetList.append(testRun)

  def runCpplint(self):
    if env['RUN_CPPLINT']:
      for path in self.cpps + self.hpps:
        lintCommand = env.Command(path + ".lint", path, lintAction)
        env.Depends(self.lib, lintCommand)
