# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import glob
import os
import subprocess
import textwrap
import SCons
from SCons.Script.SConscript import SConsEnvironment

import Helper
import SGppConfigure

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 1)
print "Using SCons", SCons.__version__

sconsVersion = SConsEnvironment()._get_major_minor_revision(SCons.__version__)
if sconsVersion < (2, 3, 0):
  Helper.printWarning("You are using an older version of scons than we do!",
                      "SGpp officially supports scons >= 2.3.0.",
                      "There are reports that it also compiles with scons >= 2.1.0.")

# to ignore folders containing a SConscript file, do the following:
# ignoreFolders = ["jsgpp"]
ignoreFolders = []

# find all modules
moduleFolders, languageSupport = Helper.getModules(ignoreFolders)

moduleNames = []
languageSupportNames = []

for moduleFolder in moduleFolders:
  moduleName = "SG_" + moduleFolder.upper()
  moduleNames.append(moduleName)

for wrapper in languageSupport:
  if wrapper == "pysgpp":
    languageSupportNames.append("SG_PYTHON")
  elif wrapper == "jsgpp":
    languageSupportNames.append("SG_JAVA")

print "Available modules:", ", ".join(moduleNames)
print "Available language support:", ", ".join(languageSupportNames)

# Define and read variables
#########################################################################

vars = Variables("custom.py")

# define the flags
vars.Add("CPPFLAGS", "Set additional compiler flags, they are compiler-dependent " +
                     "(multiple flags combined with comma, e.g. -Wall,-Wextra)", "",
                     converter=Helper.multiParamConverter)
vars.Add("LINKFLAGS", "Set additional linker flags, they are linker-dependent " +
                      "(multiple flags combined with comma, e.g. -lpython,-lm)", "",
                     converter=Helper.multiParamConverter)
vars.Add("CPPPATH", "Set path where to look for additional headers", "")
vars.Add("LIBPATH", "Set path where to look for additional libraries", "")
vars.Add("ARCH", "Set the architecture, the possible values are compiler-dependent, " +
                 "for COMPILER=gnu, e.g., the following values are possible: " +
                 "sse3, sse42, avx, fma4, avx2, avx512", "sse3")
vars.Add("COMPILER", "Set the compiler, \"gnu\" means using gcc with standard configuration, " +
                     "the following values are possible: " +
                     "gnu, clang, intel, openmpi, mpich, intel.mpi; " +
                     "when using the Intel Compiler, version 11 or higher must be used", "gnu")
vars.Add(BoolVariable("OPT", "Set compiler optimization on and off", False))
vars.Add(BoolVariable("RUN_PYTHON_TESTS", "Run Python unit tests", True))
vars.Add(BoolVariable("PYDOC", "Build Python wrapper with docstrings",
                      "SG_PYTHON" in languageSupportNames))
vars.Add(BoolVariable("SG_ALL", "Default value for the other SG_* variables; " +
                                "if True, the modules must be disabled explicitly, e.g., " +
                                "by setting SG_DATADRIVEN=0; " +
                                "if False, the modules must be enabled explicitly, e.g., " +
                                "by setting SG_DATADRIVEN=1", True))
vars.Add(BoolVariable("SG_PYTHON", "Build with Python support (default: value of SG_ALL)", None))
vars.Add(BoolVariable("SG_JAVA", "Build with Java support (default: value of SG_ALL)", None))

for moduleName in moduleNames:
  vars.Add(BoolVariable(moduleName, "Build the module " + moduleName +
                                    " (default: value of SG_ALL)", None))

vars.Add("PREFIX", "Set prefix of the paths where " +
                   "architecture-independent files are installed", "/usr/local")
vars.Add("EPREFIX", "Set prefix of the path where " +
                    "architecture-dependent files are installed (default: PREFIX)")
vars.Add("LIBDIR", "Set path where the built libraries are installed " +
                   "(default: EPREFIX/lib)")
vars.Add("INCLUDEDIR", "Set path where the header files are installed " +
                       "(default: PREFIX/include)")
vars.Add(BoolVariable("VERBOSE", "Enable verbose output", False))
vars.Add("CMD_LOGFILE", "Set path to a file to capture the build log", "build.log")
vars.Add(BoolVariable("USE_OCL", "Enable OpenCL support (only actually enabled if " +
                                 "also the OpenCL environment variables are set)", False))
vars.Add("OCL_INCLUDE_PATH", "Set path to the OpenCL header files (parent directory of CL/)")
vars.Add("OCL_LIBRARY_PATH", "Set path to the OpenCL library")
vars.Add("BOOST_INCLUDE_PATH", "Set path to the Boost header files", "/usr/include")
vars.Add("BOOST_LIBRARY_PATH", "Set path to the Boost library", "/usr/lib/x86_64-linux-gnu")
vars.Add(BoolVariable("COMPILE_BOOST_TESTS",
                      "Compile the test cases written using Boost Test", True))
vars.Add(BoolVariable("COMPILE_BOOST_PERFORMANCE_TESTS",
                      "Compile the performance tests written using Boost Test. " +
                      "Currently only buildable with OpenCL enabled", False))
vars.Add(BoolVariable("RUN_BOOST_TESTS", "Run the test cases written using Boost Test " +
                                         "(only if COMPILE_BOOST_TESTS is true)", True))
vars.Add(BoolVariable("RUN_CPPLINT",
                      "Check compliance to Google's style guide using cpplint", True))

vars.Add(BoolVariable("USE_ARMADILLO", "Set if Armadillo should be used " +
                                       "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_EIGEN", "Set if Eigen should be used " +
                                   "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_GMMPP", "Set if Gmm++ should be used " +
                                   "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_UMFPACK", "Set if UMFPACK should be used " +
                                     "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("BUILD_STATICLIB", "Set if static libraries should be built " +
                                         "instead of shared libraries", False))
vars.Add(BoolVariable("PRINT_INSTRUCTIONS", "Print instructions for installing SG++", True))

# create temporary environment to check which system and compiler we should use
# (the Environment call without "tools=[]" crashes with MinGW,
# so we do it like that)
env = Environment(variables=vars, ENV=os.environ, tools=[])

if (env["PLATFORM"].lower() == "win32") and \
   (env["COMPILER"].lower() == "gnu"):
  # MinGW: use gcc toolschain
  tools = ["gnulink", "gcc", "g++", "gas", "ar", "swig"]
else:
  # otherwise: use default toolchain
  tools = ["default"]

# Initialize environment
#########################################################################

env = Environment(variables=vars, ENV=os.environ, tools=tools)

# fail if unknown variables where encountered on the command line
unknownVariables = [var for var in vars.UnknownVariables()
                    if var not in ["CXX", "CC", "CFLAGS", "CPPDEFINES"]]
if len(unknownVariables) > 0:
  Helper.printErrorAndExit("The following command line variables could not be recognized:",
                           unknownVariables,
                           "Type \"scons -h\" to list all available variables.")

# set "dynamic" defaults of variables
# (which depend on the values of other variables)
env["SG_PYTHON"] = env.get("SG_PYTHON",
                           env["SG_ALL"] and ("SG_PYTHON" in languageSupportNames))
env["SG_JAVA"] = env.get("SG_JAVA",
                         env["SG_ALL"] and ("SG_JAVA" in languageSupportNames))

for moduleName in moduleNames:
  env[moduleName] = env.get(moduleName, env["SG_ALL"])

env["EPREFIX"] = env.get("EPREFIX", env["PREFIX"])
env["LIBDIR"] = env.get("LIBDIR", os.path.join(env["EPREFIX"], "lib"))
env["INCLUDEDIR"] = env.get("INCLUDEDIR", os.path.join(env["PREFIX"], "include"))

# don't create the Doxyfile if we're cleaning:
if not env.GetOption("clean"):
  Helper.prepareDoxyfile(moduleFolders)

if "CXX" in ARGUMENTS:
  print "CXX: ", ARGUMENTS["CXX"]
  env["CXX"] = ARGUMENTS["CXX"]
if "CC" in ARGUMENTS:
  print "CC: ", ARGUMENTS["CC"]
  env["CC"] = ARGUMENTS["CC"]
if "CPPFLAGS" in ARGUMENTS:
  env["CPPFLAGS"] = ARGUMENTS["CPPFLAGS"].split(",")
if "CFLAGS" in ARGUMENTS:
  env["CFLAGS"] = ARGUMENTS["CFLAGS"]
env.AppendUnique(CPPDEFINES = {})
if "CPPDEFINES" in ARGUMENTS:
  for define in ARGUMENTS["CPPDEFINES"].split(","):
    key, value = define.split("=")
    env["CPPDEFINES"][key] = value

if "CPPPATH" in ARGUMENTS:
  env["CPPPATH"] = ARGUMENTS["CPPPATH"].split(",")
if "LIBPATH" in ARGUMENTS:
  env["LIBPATH"] = ARGUMENTS["LIBPATH"].split(",")

env.Export("moduleNames")
env.Export("moduleFolders")

# help text
Help("""
--------------------------------------------------------------------------------

Type "scons" to build the libraries.

You can set build options in the command line. For example, use
"scons COMPILER=intel" to use the Intel compiler. Boolean flags such as OPT can
be set using various notations (1/0, true/false, yes/no). Use the SG_* flags to
disable building of specific modules (e.g., "scons SG_DATADRIVEN=0"). If you
want to see the compiler and linker commands, use "scons VERBOSE=1". Use
"scons -j 4" to enable parallel building on four cores.

You can find a list of options below this text.

Type "scons doxygen" to build the Doxygen documentation, which can afterwards be
accessed via doc/html/index.html. If you additionally exclude modules (e.g.,
"scons SG_DATADRIVEN=0 doxygen"), then only the documentation for the remaining
modules will be built.

--------------------------------------------------------------------------------

Parameters are:
""" +
"\n".join([textwrap.fill(line, 80)
           for line in vars.GenerateHelpText(env).splitlines()]))

# add trailing slashes were required and if not present
BUILD_DIR = Dir(os.path.join("lib", "sgpp"))
Export("BUILD_DIR")
PYSGPP_PACKAGE_PATH = Dir(os.path.join("lib"))
Export("PYSGPP_PACKAGE_PATH")
PYSGPP_BUILD_PATH = Dir(os.path.join(PYSGPP_PACKAGE_PATH.abspath, "pysgpp"))
Export("PYSGPP_BUILD_PATH")
JSGPP_BUILD_PATH = Dir(os.path.join("lib", "jsgpp"))
Export("JSGPP_BUILD_PATH")
EXAMPLE_DIR = Dir(os.path.join("bin", "examples"))
Export("EXAMPLE_DIR")

# don't configure if we're cleaning:
if not env.GetOption("clean"):
  SGppConfigure.doConfigure(env, moduleFolders, languageSupport)

# fix for "command line too long" errors on MinGW
# (from https://bitbucket.org/scons/scons/wiki/LongCmdLinesOnWin32)
if env["PLATFORM"] == "win32":
  Helper.setWin32Spawn(env)

# add #/lib/sgpp to LIBPATH
# (to add corresponding -L... flags to linker calls)
env.Append(LIBPATH=[BUILD_DIR])

# add C++ defines for all modules
for module in moduleNames:
  if env[module]:
    env["CPPDEFINES"][module] = "1"

# environement setup finished, export environment
Export("env")

env.Append(CPPPATH=["#/tools"])
config = env.Configure()
Export("config")

# update PATH under win32/LD_LIBRARY_PATH otherwise
# to add BUILD_DIR (so we can run the Boost tests)
if env["PLATFORM"] == "win32":
  env["ENV"]["PATH"] = os.pathsep.join([env["ENV"].get("PATH", ""),
                                        BUILD_DIR.abspath])

  # also add the Boost library path to the PATH
  # so that the Boost test *.dll can be found when running the tests
  if env["RUN_BOOST_TESTS"]:
    env["ENV"]["PATH"] = os.pathsep.join([env["ENV"].get("PATH", ""),
                                          env["BOOST_LIBRARY_PATH"]])
# Mac OS X doesn't use LD_LIBRARY_PATH
elif env["PLATFORM"] == "darwin":
  env["ENV"]["DYLD_FALLBACK_LIBRARY_PATH"] = os.pathsep.join([
      env["ENV"].get("DYLD_FALLBACK_LIBRARY_PATH", ""),
      BUILD_DIR.abspath])
else:
  env["ENV"]["LD_LIBRARY_PATH"] = os.pathsep.join([
      env["ENV"].get("LD_LIBRARY_PATH", ""),
      BUILD_DIR.abspath])

# Add the pysgpp package path to the environment
#########################################################################

if env["PLATFORM"] == "win32":
  # try to import pysgpp to detect an already existing installation, which
  # could cause trouble
  try:
    import pysgpp
    Helper.printWarning("An existing installation of pysgpp was detected."
                        "To get rid of this warning remove the pysgpp package"
                        "from your local Oython installation.")
  except:
    pass

  # get a temporary folders
  import tempfile, uuid
  # get temp directory
  pysgppTempFolder = os.path.join(tempfile.gettempdir(),
                                  "site-pyspp-" + str(uuid.uuid4()))
  # create temp folder
  os.makedirs(pysgppTempFolder)

  # add it to the build python path
  env["ENV"]["PYTHONPATH"] = os.pathsep.join([pysgppTempFolder,
                                              env["ENV"].get("PYTHONPATH", "")])
else:
  env["ENV"]["PYTHONPATH"] = os.pathsep.join([env["ENV"].get("PYTHONPATH", ""),
                                              PYSGPP_PACKAGE_PATH.abspath])

# Style checker
#########################################################################

def lintAction(target, source, env):
  p = subprocess.Popen(["python", "tools/cpplint.py", "--ignorecfg=yes",
                        "--extensions=cpp,hpp", "--linelength=100",
                        source[0].abspath],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  # wait for termination and get output on stdout and stderr
  stdout, stderr = p.communicate()
  # cpplint prints on stderr
  for line in stderr.splitlines():
    # skip status lines, empty lines, and some warning types
    if ("Done processing " in line) or \
        ("Total errors found: " in line) or \
        ("Is this a non-const reference? " +
         "If so, make const or use a pointer:" in line) or \
        ("Consider using rand_r(...) instead of rand(...) for " +
         "improved thread safety." in line) or \
        ("<chrono> is an unapproved C++11 header." in line) or \
        (line == ""):
      pass
    else:
      parts = line.split(":  ")
      location = parts[0]
      message = ":  ".join(parts[1:])
      print location + ": warning: " + message
  # touch file without writing anything
  # (to indicate for the next run of SCons that we already checked this file)
  with open(target[0].abspath, "w"): pass

env.Export("lintAction")

# Custom builders for Python and Boost tests
#########################################################################

if env["RUN_PYTHON_TESTS"] and env["SG_PYTHON"]:
  # do the actual thing
  builder = Builder(action="python $SOURCE", chdir=0)
  env.Append(BUILDERS={"Test" : builder})
  builder = Builder(action="python $SOURCE")
  env.Append(BUILDERS={"SimpleTest" : builder})

if env["COMPILE_BOOST_TESTS"]:
  builder = Builder(action="./$SOURCE")
  env.Append(BUILDERS={"BoostTest" : builder})

# Building the modules
#########################################################################

libraryTargetList = []
pythonTestTargetList = []
boostTestTargetList = []
boostTestRunTargetList = []
exampleTargetList = []
pydocTargetList = []
headerSourceList = []
headerDestList = []
env.Export("libraryTargetList")
env.Export("pythonTestTargetList")
env.Export("boostTestTargetList")
env.Export("boostTestRunTargetList")
env.Export("exampleTargetList")
env.Export("pydocTargetList")
env.Export("headerSourceList")
env.Export("headerDestList")

# compile selected modules
flattenedDependencyGraph = []

for moduleFolder in moduleFolders:
  if not env["SG_" + moduleFolder.upper()]:
    continue

  print "Preparing to build module:", moduleFolder
  # SConscript("src/sgpp/SConscript" + moduleFolder, variant_dir="#/tmp/build/", duplicate=0)
  env.SConscript("#/" + moduleFolder + "/SConscript", {"env": env, "moduleName": moduleFolder})

  # add the dependencies of the current module to the overall dependency graph
  Import("moduleDependencies")
  Import("libname")
  flattenedDependencyGraph = Helper.flatDependencyGraph(
      [libname] + moduleDependencies, flattenedDependencyGraph)

Export("flattenedDependencyGraph")

# compile pysgpp
if env["SG_PYTHON"]:
  env.SConscript("#/pysgpp/SConscript", {"env": env, "moduleName": "pysgpp"})

# compile jsgpp
if env["SG_JAVA"]:
  env.SConscript("#/jsgpp/SConscript", {"env": env, "moduleName": "jsgpp"})

if (not env["SG_PYTHON"]) and (not env["SG_JAVA"]) and (len(flattenedDependencyGraph) == 0):
  Helper.printErrorAndExit("You must enable at least one module (e.g., SG_BASE=1).")

# Python tests
#########################################################################

finalStepDependencies = [libraryTargetList, pydocTargetList]

# prepare python package for unit testing
def installPythonLibToTmp(target, source, env):
  # get temp directory
  pysgppTempFolder = source[0].get_string(0)

  # install python interface to tmp directory
  p = subprocess.call(["python", "setup.py",
                       "--quiet",
                       "install", "--install-lib=%s" % pysgppTempFolder])
  if p != 0:
    Helper.printErrorAndExit("Installing Python package to the temporary folder",
                             pysgppTempFolder,
                             "failed. Python unit tests cannot be run automatically.")

if env["RUN_PYTHON_TESTS"] and env["SG_PYTHON"]:
  if env["PLATFORM"] == "win32":
    # install the python library to that temporary folder
    installPythonLibToTmpCommand = env.Command("installPythonLibToTmp", [pysgppTempFolder],
                                               installPythonLibToTmp)
    env.Depends(installPythonLibToTmpCommand, finalStepDependencies)
    env.Depends(pythonTestTargetList, installPythonLibToTmpCommand)

  # serialize tests (from https://bitbucket.org/scons/scons/wiki/FrequentlyAskedQuestions#
  # markdown-header-how-do-i-prevent-commands-from-being-executed-in-parallel)
  # and move them at the end of the build
  env.Depends(pythonTestTargetList, finalStepDependencies)
  finalStepDependencies.append(pythonTestTargetList)
  env.SideEffect("sideEffectFinalSteps", pythonTestTargetList)

# Boost tests
#########################################################################

if env["COMPILE_BOOST_TESTS"]:
  env.Depends(boostTestTargetList, finalStepDependencies)
  finalStepDependencies.append(boostTestTargetList)
  env.SideEffect("sideEffectFinalSteps", boostTestTargetList)

  if env["RUN_BOOST_TESTS"]:
    env.Depends(boostTestRunTargetList, finalStepDependencies)
    finalStepDependencies.append(boostTestRunTargetList)
    env.SideEffect("sideEffectFinalSteps", boostTestRunTargetList)

# Examples
#########################################################################

env.Depends(exampleTargetList, finalStepDependencies)
finalStepDependencies.append(exampleTargetList)
env.SideEffect("sideEffectFinalSteps", exampleTargetList)

# Final output
#########################################################################

def printInstructions(target, source, env):
  import string
  if env["PLATFORM"] in ["cygwin", "win32"]:
    filename = "INSTRUCTIONS_WINDOWS"
  elif env["PLATFORM"] == "darwin" :
    filename = "INSTRUCTIONS_MAC"
  else:
    filename = "INSTRUCTIONS"

  with open(filename) as f:
    instructionsTemplate = string.Template(f.read())
    print
    print instructionsTemplate.safe_substitute(SGPP_BUILD_PATH=BUILD_DIR.abspath,
                                               PYSGPP_PACKAGE_PATH=PYSGPP_PACKAGE_PATH.abspath)

printInstructionsTarget = env.Command("printInstructions", [], printInstructions)
env.Depends(printInstructionsTarget, finalStepDependencies)
finalStepDependencies.append(printInstructionsTarget)

# System-wide installation
#########################################################################

installLibSGpp = env.Alias("install-lib-sgpp",
                           env.Install(os.path.join(env.get("LIBDIR"), "sgpp"),
                                       libraryTargetList))

headerFinalDestList = []
for headerDest in headerDestList:
  headerFinalDestList.append(os.path.join(env.get("INCLUDEDIR"), headerDest))

installIncSGpp = env.Alias("install-inc-sgpp",
                           env.InstallAs(headerFinalDestList, headerSourceList))

env.Alias("install", [installLibSGpp, installIncSGpp])

# Doxygen
#########################################################################

doxygen = env.Command("doc/xml/index.xml", "Doxyfile", "doxygen $SOURCE")
env.Alias("doxygen", doxygen)

# Things to be cleaned
#########################################################################

# build + log files
env.Clean("clean", ["lib/", "jsgpp/java/", "config.log", "build.log"])
# Doxygen stuff
env.Clean("clean", ["Doxyfile", "doxygen_warnings.log", "doc/html/", "doc/xml/",
                    "base/doc/doxygen/examples.doxy", "base/doc/doxygen/modules.doxy"])

for module in moduleFolders:
  # PYDOC stuff
  env.Clean("clean", [module + "/Doxyfile", module + "/doc/doxygen_warnings.log",
                      module + "/doc/xml/"])

# Default targets
#########################################################################

if not GetOption("clean"):
  env.Default(finalStepDependencies)
else:
  env.Default(finalStepDependencies + ["clean"])
