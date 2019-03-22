# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import fnmatch
import inspect
import os

import Helper

class Module(object):
  def __init__(self, moduleDependencies, additionalDependencies=[],
               additionalBoostTestDependencies=[], excludeFiles=[]):
    # moduleDependencies: list of SG++ modules (e.g., "sgppbase") which are needed for this module
    self.moduleDependencies = moduleDependencies
    # additionalDependencies: list of other libraries on which the module depends
    # (used for the module library itself, the examples, and the Boost tests)
    self.additionalDependencies = additionalDependencies
    # additionalBoostTestDependencies: list of other libraries on which the Boost tests depend
    self.additionalBoostTestDependencies = additionalBoostTestDependencies
    # excludeFiles: source/header files which should be excluded from building
    self.excludeFiles = excludeFiles
    self.cpps = []
    self.hpps = []
    self.objs = []

    # inject those variables in the module namespace
    # which were imported in the SConscript files, including env;
    # this assumes that the __init__ was called directly in the SConscript at top-level,
    # not in a function
    # (not a very nice method, but passing the variables in function arguments
    # would be very tedious)
    variables = inspect.stack()[1][0].f_globals
    for name, value in variables.items():
      globals()[name] = value

    # check module dependencies, fail if a dependency is missing
    for module in self.moduleDependencies:
      if module.startswith("sgpp") and (not env["SG_" + module[4:].upper()]):
        Helper.printErrorAndExit(
            "The \"{}\" module depends on the \"{}\" module, ".format(moduleName, module[4:]),
            "which is currently not enabled. Please enable it by setting SG_{}=1.".format(
                module[4:].upper()))

  def scanSource(self, sourceFolder="src"):
    """Scan the given directory for source and header files.
    """
    # the sourceFolder is relative to the current directory
    # ==> make absolute path
    sourceFolder = os.path.join(Dir(".").abspath, sourceFolder)

    for currentFolder, subdirNames, fileNames in os.walk(sourceFolder, topdown=True):
      # sort subdirectories and files
      subdirNames.sort()
      fileNames.sort()

      # check if there's a submodule SConscript in the current folder
      sconscriptPath = os.path.join(currentFolder, "SConscript")

      if (currentFolder != sourceFolder) and os.path.exists(sconscriptPath):
        # call submodule SConscript, tell it that's it's our module calling
        module = self
        env.SConscript(sconscriptPath, exports="module")
        # remove subfolders from iteration as they are already processed
        # (this is why topdown=True is required)
        subdirNames[:] = []
      else:
        # process source files
        for fileName in fnmatch.filter(fileNames, "*.cpp"):
          if fileName in self.excludeFiles: continue
          cpp = os.path.join(currentFolder, fileName)
          self.cpps.append(cpp)
          self.objs.append(env.SharedObject(cpp))

          # process cuda source files if cuda is enabled
        if env["USE_CUDA"]:
          for fileName in fnmatch.filter(fileNames, "*.cu"):
            if fileName in self.excludeFiles: continue
            cpp = os.path.join(currentFolder, fileName)
            self.cpps.append(cpp)
            self.objs.append(env.SharedObject(cpp))

        # process header files
        for fileName in fnmatch.filter(fileNames, "*.hpp"):
          if fileName in self.excludeFiles: continue
          hpp = os.path.join(currentFolder, fileName)
          self.hpps.append(hpp)

    # append headers to install list
    for hpp in self.hpps:
      headerSourceList.append(hpp)
      headerDestList.append(hpp)

  def buildLibrary(self):
    """Build the module.
    """

    # name of the library to build
    self.libname = "sgpp" + moduleName

    # change library names if we're linking statically
    if env["BUILD_STATICLIB"]:
      self.libname += "static"
      self.moduleDependencies = [module + "static" for module in self.moduleDependencies]

    # export library name and dependencies
    libname = self.libname
    env.Export("libname")
    moduleDependencies = self.moduleDependencies
    env.Export("moduleDependencies")

    if env["BUILD_STATICLIB"]:
      # build static library
      libsuffix = env["LIBSUFFIX"]
      envClone = env.Clone()
      envClone.AppendUnique(LIBS=self.moduleDependencies + self.additionalDependencies)
      self.lib = envClone.StaticLibrary(target=self.libname, source=self.objs)
    else:
      # build shared library
      libsuffix = env["SHLIBSUFFIX"]
      envClone = env.Clone()
      envClone.AppendUnique(LIBS=self.moduleDependencies + self.additionalDependencies)
      self.lib = envClone.SharedLibrary(target=self.libname, source=self.objs)

    # set module dependencies
    for module in self.moduleDependencies:
      if module.startswith("sgpp"):
        otherLib = os.path.join("#", BUILD_DIR.path, env["LIBPREFIX"] + module + libsuffix)
        env.Depends(self.lib, otherLib)

    # install the library
    self.libInstall = env.Install(BUILD_DIR, self.lib)
    libraryTargetList.append(self.libInstall)

  def generatePythonDocstrings(self):
    """Generate Python docstrings for pysgpp.
    """
    if env["PYDOC"] and env["SG_PYTHON"]:
      # read module Doxygen file
      with open(File("#/moduleDoxy").abspath, "r") as f:
        moduleDoxy = f.read()

      if not env.GetOption("clean"):
        # write module-specific Doxyfile
        with open("Doxyfile", "w") as f:
          f.write(moduleDoxy.replace("$modname", moduleName).replace("$quiet", "YES"))

      # execute Doxygen
      doxygen = env.Command("doc/xml/index.xml", "Doxyfile", "doxygen $SOURCE")
      env.Depends(doxygen, self.cpps + self.hpps)

      # call doxy2swig.py, converting the Doxygen XML output to a SWIG *.i file
      doxy2swig = env.Command(os.path.join("#", "pysgpp", moduleName + "_doc.i"), doxygen,
                              "python pysgpp/doxy2swig.py -o -c -q $SOURCE $TARGET")
      pydocTargetList.append(doxy2swig)

  def buildExamples(self, exampleFolder="examples", additionalExampleDependencies=[]):
    """Build the examples.
    """
    # set libraries
    exampleEnv = env.Clone()
    exampleEnv.AppendUnique(LIBS=[self.libname] +
                                 self.moduleDependencies + self.additionalDependencies +
                                 additionalExampleDependencies)

    # for each example
    for fileName in os.listdir(exampleFolder):
      if fnmatch.fnmatch(fileName, "*.cpp"):
        # source file
        cpp = os.path.join(exampleFolder, fileName)
        self.cpps.append(cpp)
        example = exampleEnv.Program(source=cpp)
        exampleEnv.Depends(example, self.libInstall)
        exampleTargetList.append(example)
      elif fnmatch.fnmatch(fileName, "*.hpp"):
        # header file
        hpp = os.path.join(exampleFolder, fileName)
        self.hpps.append(hpp)

  def runExamples(self, exampleFolder="examples", language="all"):
    """Run the examples.
    """
    if language == "all":
      for language in ["cpp", "python"]:
        self.runExamples(exampleFolder="examples", language=language)
      return
    elif language == "cpp":
      if not env["RUN_CPP_EXAMPLES"]: return
      fileNameFilter = "*.cpp"
    elif language == "python":
      if not env["RUN_PYTHON_EXAMPLES"]: return
      fileNameFilter = "*.py"
    else:
      raise ValueError("Unsupported language for running examples.")

    for fileName in os.listdir(exampleFolder):
      if fnmatch.fnmatch(fileName, fileNameFilter):
        sourcePath = os.path.join(exampleFolder, fileName)

        if language == "cpp":
          sourcePath = os.path.splitext(sourcePath)[0]
          cppExampleRun = env.CppExample(sourcePath + "_run", sourcePath)
          cppTestRunTargetList.append(cppExampleRun)
        elif language == "python":
          pythonExampleRun = env.PythonExample(sourcePath + "_run", sourcePath)
          pythonTestRunTargetList.append(pythonExampleRun)

  def runPythonTests(self):
    """Run the Python tests.
    """
    if env["RUN_PYTHON_TESTS"] and env["SG_PYTHON"]:
      # run Python test
      moduleTest = env.Test(os.path.join("tests", "test_{}.py".format(moduleName)))
      pythonTestTargetList.append(moduleTest)

  def buildBoostTests(self, boostTestFolder="tests", compileFlag="COMPILE_BOOST_TESTS"):
    """Compile the Boost tests.
    """
    if env[compileFlag]:
      # set libraries
      testEnv = env.Clone()
      testEnv.AppendUnique(LIBS=[self.libname] +
                                self.moduleDependencies + self.additionalDependencies +
                                ["boost_unit_test_framework"] +
                                self.additionalBoostTestDependencies)

      testObjs = []

      # for each test
      for currentFolder, subdirNames, fileNames in os.walk(boostTestFolder, topdown=True):
        for fileName in fnmatch.filter(fileNames, "*.cpp"):
          # source file
          cpp = os.path.join(currentFolder, fileName)
          self.cpps.append(cpp)
          testObjs.append(testEnv.SharedObject(cpp))
        for fileName in fnmatch.filter(fileNames, "*.hpp"):
          # header file
          hpp = os.path.join(currentFolder, fileName)
          self.hpps.append(hpp)

      # only build Boost test executable if there are any tests
      if len(testObjs) > 0:
        self.boostTestExecutable = \
            os.path.join(boostTestFolder, "test_{}_boost".format(moduleName)) + \
            (".exe" if env["PLATFORM"] == "win32" else "")
        test = testEnv.Program(self.boostTestExecutable, testObjs)
        testEnv.Depends(test, self.libInstall)
        boostTestTargetList.append(test)

  def runBoostTests(self, boostTestFolder="tests",
                    compileFlag="COMPILE_BOOST_TESTS", runFlag="RUN_BOOST_TESTS"):
    """Run the Boost tests.
    """
    if env[compileFlag] and env[runFlag]:
      # run Boost tests
      testRun = env.BoostTest(self.boostTestExecutable + "_run", self.boostTestExecutable)
      boostTestRunTargetList.append(testRun)

  def runCpplint(self):
    """Run the style checker.
    """
    if env['RUN_CPPLINT']:
      # run the style checker on all source and header files
      for path in self.cpps + self.hpps:
        lintCommand = env.Command(path + ".lint", path, lintAction)
        env.Depends(self.lib, lintCommand)
