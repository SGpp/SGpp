#! /usr/bin/env python3
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import atexit
import glob
import os
import pipes
import platform
import subprocess
import sys
import textwrap
import SCons
from SCons.Script.SConscript import SConsEnvironment

import DoxygenHelper
import Helper
import SGppConfigure

# detour stdout and stderr to file
sys.stdout = Helper.Logger(sys.stdout)
sys.stderr = Helper.Logger(sys.stderr)

# object responsible for printing a final message depending on success or error
# with additional instructions
finalMessagePrinter = Helper.FinalMessagePrinter()
atexit.register(finalMessagePrinter.printMessage)

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 1)
Helper.printInfo("Platform: {}".format(", ".join(platform.uname())))
Helper.printInfo("Using SCons {} on Python {}".format(SCons.__version__, platform.python_version()))
Helper.printInfo("SCons command line: {}".format(" ".join([pipes.quote(arg) for arg in sys.argv])))

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
  elif wrapper == "matsgpp":
    languageSupportNames.append("SG_MATLAB")

Helper.printInfo("Available modules: {}".format(", ".join(moduleNames)))
Helper.printInfo("Available language support: {}".format(", ".join(languageSupportNames)))

# Define and read variables
#########################################################################

AddOption('--configfile',
          dest='configfile',
          type='string',
          nargs=1,
          action='store',
          help='specify custom options input file (see custom.py in scons documentation)')
input_file = GetOption('configfile')

if input_file:
  vars = Variables(input_file)
else:
  vars = Variables("custom.py")


# define the flags
vars.Add("CFLAGS", "Set additional C compiler flags, they are compiler-dependent "
                   "(multiple flags separated by space: '-Wall -Wextra')", "",
                   converter=Helper.multiParamConverter)
vars.Add("CPPFLAGS", "Set additional C++ compiler flags, they are compiler-dependent "
                     "(multiple flags separated by space: '-Wall -Wextra')", "",
                     converter=Helper.multiParamConverter)
vars.Add("LINKFLAGS", "Set additional linker flags, they are linker-dependent "
                      "(multiple flags separated by space: '-lpython -lm')", "",
                     converter=Helper.multiParamConverter)
vars.Add("CPPDEFINES", "Set additional C++ defines "
                       "(multiple defines separated by space: 'FLAG_A=1 FLAG_B=2')", "",
                       converter=Helper.multiParamDefineConverter)
vars.Add("CPPPATH", "Set path where to look for additional headers "
                    "(multiple paths separated by '{}': '{}')".format(
                        os.pathsep, os.pathsep.join([os.path.join("first", "path"),
                                                     os.path.join("second", "path")])), "",
                    converter=Helper.multiParamPathConverter)
vars.Add("LIBPATH", "Set path where to look for additional libraries", "",
                    converter=Helper.multiParamPathConverter)
vars.Add("ARCH", "Set the architecture, the possible values are compiler-dependent, " +
                 "for COMPILER=gnu, e.g., the following values are possible: " +
                 "sse3, sse42, avx, fma4, avx2, avx512", "sse3")
vars.Add("COMPILER", "Set the compiler, \"gnu\" means using gcc with standard configuration, " +
                     "the following values are possible: " +
                     "gnu, clang, intel, openmpi, mpich, intel.mpi; " +
                     "when using the Intel Compiler, version 11 or higher must be used", "gnu")
vars.Add("CC", "Override the C compiler, can be used to select a specific compiler version, otherwise use \"COMPILER\"", None)
vars.Add("CXX", "Override the C++ compiler, can be used to select a specific compiler version, otherwise use \"COMPILER\"", None)
vars.Add("LINK", "Override the linker, can be used to select a specific linker version", None)
vars.Add(BoolVariable("OPT", "Set compiler optimization on and off", True))
vars.Add(BoolVariable("RUN_ON_HAZELHEN", "Add some special options on hazelhen", False))
vars.Add(BoolVariable("RUN_PYTHON_TESTS", "Run Python unit tests", True))
vars.Add(BoolVariable("PYDOC", "Build Python wrapper with docstrings", False))
vars.Add(BoolVariable("USE_PYTHON2_FOR_PYSGPP", "Enforce using Python 2.x for pysgpp", False))
vars.Add(BoolVariable("SG_ALL", "Default value for the other SG_* variables; " +
                                "if True, the modules must be disabled explicitly, e.g., " +
                                "by setting SG_DATADRIVEN=0; " +
                                "if False, the modules must be enabled explicitly, e.g., " +
                                "by setting SG_DATADRIVEN=1", True))
vars.Add(BoolVariable("SG_PYTHON", "Build with Python support (default: value of SG_ALL)", None))
vars.Add(BoolVariable("SG_JAVA", "Build with Java support (default: value of SG_ALL)", None))
vars.Add(BoolVariable("SG_MATLAB", "Build with MATLAB support", False))
vars.Add("SWIGFLAGS", "Set additional SWIG flags, they are compiler-dependent "
                      "(multiple flags separated by space: '-Wall -Wextra')", "",
                      converter=Helper.multiParamConverter)

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
vars.Add(BoolVariable("VERBOSE", "Enable verbose output", True))
vars.Add(BoolVariable("USE_OCL", "Enable OpenCL support (only actually enabled if " +
                                 "also the OpenCL environment variables are set)", False))
vars.Add(BoolVariable("USE_CUDA", "Enable CUDA support (you might need to provide an 'CUDA_TOOLKIT_PATH')", False))
vars.Add(BoolVariable("USE_HPX", "Enable HPX support (implies USE_OCL)", False))
vars.Add("OCL_INCLUDE_PATH", "Set path to the OpenCL header files (parent directory of CL/)")
vars.Add("OCL_LIBRARY_PATH", "Set path to the OpenCL library")
vars.Add("BOOST_INCLUDE_PATH", "Set path to the Boost header files", "/usr/include")
vars.Add("BOOST_LIBRARY_PATH", "Set path to the Boost library", None)
vars.Add("MATLAB_INCLUDE_PATH", "Set path to the MATLAB header files", None)
vars.Add("MATLAB_LIBRARY_PATH", "Set path to the MATLAB libraries", None)
vars.Add("HPX_DEBUG_LIBRARY_PATH", "Sets the path to the HPX debug libraries", None)
vars.Add("HPX_RELEASE_LIBRARY_PATH", "Sets the path to the HPX release libraries", None)
vars.Add("HPX_SHARED_INCLUDE_PATH", "Sets the path to the HPX shared headers", None)
vars.Add("HPX_DEBUG_INCLUDE_PATH", "Sets the path to the HPX debug headers", None)
vars.Add("HPX_RELEASE_INCLUDE_PATH", "Sets the path to the HPX release headers", None)
vars.Add("GSL_INCLUDE_PATH", "Set path to the GSL header files", "/usr/include")
vars.Add("GSL_LIBRARY_PATH", "Set path to the GSL library", None)
vars.Add("SCALAPACK_LIBRARY_PATH", "Set the path to the ScaLAPACK/MKL library", None)
vars.Add("SCALAPACK_LIBRARY_NAME", "Set the name of the ScaLAPACK library", None)
vars.Add(BoolVariable("COMPILE_BOOST_TESTS",
                      "Compile the test cases written using Boost Test", True))
vars.Add(BoolVariable("COMPILE_BOOST_PERFORMANCE_TESTS",
                      "Compile the performance tests written using Boost Test. " +
                      "Currently only buildable with OpenCL enabled", False))
vars.Add(BoolVariable("RUN_BOOST_PERFORMANCE_TESTS", "Run the test cases written using Boost Test " +
                                         "(only if COMPILE_BOOST_PERFORMANCE_TESTS is true)", True))
vars.Add(BoolVariable("RUN_BOOST_TESTS", "Run the test cases written using Boost Test " +
                                         "(only if COMPILE_BOOST_TESTS is true)", True))
vars.Add(BoolVariable("CHECK_STYLE",
                      "Check compliance to Google's style guide using cpplint", True))
vars.Add(BoolVariable("RUN_CPP_EXAMPLES", "Run all C++ examples", False))
vars.Add(BoolVariable("RUN_PYTHON_EXAMPLES", "Run all Python examples", False))

vars.Add(BoolVariable("USE_ARMADILLO", "Set if Armadillo should be used " +
                                       "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_EIGEN", "Set if Eigen should be used " +
                                   "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_GMMPP", "Set if Gmm++ should be used " +
                                   "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_UMFPACK", "Set if UMFPACK should be used " +
                                     "(only relevant for sgpp::optimization)", False))
vars.Add(BoolVariable("USE_DAKOTA", "Set if Dakota library should be used " +
                                   "(only relevant for sgpp::combigrid)", False))
vars.Add(BoolVariable("USE_GSL", "Set if GNU Scientific Library should be used " +
                                     "(only relevant for sgpp::datadriven)", False))
vars.Add(BoolVariable("USE_CGAL", "Set if Computational Geometry Algorithms Library should be used " +
                                     "(only relevant for new_sgde)", False))

vars.Add(BoolVariable("USE_ZLIB", "Set if zlib should be used " +
                                     "(relevant for sgpp::datadriven to read compressed dataset files), not available for windows", False))
vars.Add(BoolVariable("USE_SCALAPACK", "Set if the ScaLAPACK library should be used " +
                                          "(requires MPI, only relevant for sgpp::datadriven)", None))
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
finalMessagePrinter.env = env

# USE_HPX implies USE_OCL
if env["USE_HPX"]:
  env["USE_OCL"] = True

# fail if unknown variables where encountered on the command line
unknownVariables = vars.UnknownVariables()
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
env["SG_MATLAB"] = env.get("SG_MATLAB",
                           env["SG_ALL"] and ("SG_MATLAB" in languageSupportNames))

for moduleName in moduleNames:
  env[moduleName] = env.get(moduleName, env["SG_ALL"])

env["EPREFIX"] = env.get("EPREFIX", env["PREFIX"])
env["LIBDIR"] = env.get("LIBDIR", os.path.join(env["EPREFIX"], "lib"))
env["INCLUDEDIR"] = env.get("INCLUDEDIR", os.path.join(env["PREFIX"], "include"))
env["BOOST_LIBRARY_PATH"] = env.get("BOOST_LIBRARY_PATH", "/usr/lib/x86_64-linux-gnu"
                                    if env["PLATFORM"] not in ["darwin", "win32"]
                                    else "")

# only create the Doxyfile if building Doxygen:
if ("doxygen" in BUILD_TARGETS) and (not env.GetOption("clean")):
  Helper.printInfo("Building Doxyfile for modules: "+
                   ', '.join([moduleFolder for moduleFolder in moduleFolders if env["SG_" + moduleFolder.upper()]]))
  DoxygenHelper.prepareDoxygen([moduleFolder for moduleFolder in moduleFolders if env["SG_" + moduleFolder.upper()]])

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
BUILD_DIR = Dir(os.path.join("lib"))
Export("BUILD_DIR")
PYSGPP_PACKAGE_PATH = Dir(os.path.join("lib"))
Export("PYSGPP_PACKAGE_PATH")
PYSGPP_BUILD_PATH = Dir(os.path.join(PYSGPP_PACKAGE_PATH.abspath, "pysgpp"))
Export("PYSGPP_BUILD_PATH")
JSGPP_BUILD_PATH = Dir(os.path.join("lib", "jsgpp"))
Export("JSGPP_BUILD_PATH")
MATSGPP_BUILD_PATH = Dir(os.path.join("lib", "matsgpp"))
Export("MATSGPP_BUILD_PATH")
EXAMPLE_DIR = Dir(os.path.join("bin", "examples"))
Export("EXAMPLE_DIR")

# save ARGUMENTS dict (unparsed command-line arguments) in env
# in order to be able to tell, e.g., if the user has specified a custom value
# for CXX or if the default ("g++" on Linux) has been used
env.arguments = ARGUMENTS

if not ( env.GetOption('clean') or env.GetOption('help') ):
  SGppConfigure.doConfigure(env, moduleFolders, languageSupport)

# fix for "command line too long" errors on MinGW
# (from https://bitbucket.org/scons/scons/wiki/LongCmdLinesOnWin32)
if env["PLATFORM"] == "win32":
  Helper.setWin32Spawn(env)
else:
  Helper.setSpawn(env)

# add #/lib to LIBPATH
# (to add corresponding -L... flags to linker calls)
env.Append(LIBPATH=[BUILD_DIR])

# environement setup finished, export environment
Export("env")

env.Append(CPPPATH=["#/tools"])
config = env.Configure()
Export("config")

# update PATH under win32/LD_LIBRARY_PATH otherwise
# to add BUILD_DIR (so we can run the Boost tests)
if env["PLATFORM"] == "win32":
  env["ENV"]["PATH"] = os.pathsep.join([BUILD_DIR.abspath,
                                        env["ENV"].get("PATH", "")])

  # also add the Boost library path to the PATH
  # so that the Boost test *.dll can be found when running the tests
  if env["COMPILE_BOOST_TESTS"]:
    env["ENV"]["PATH"] = os.pathsep.join([env["BOOST_LIBRARY_PATH"],
                                          env["ENV"].get("PATH", "")])
# Mac OS X doesn't use LD_LIBRARY_PATH
elif env["PLATFORM"] == "darwin":
  env["ENV"]["DYLD_FALLBACK_LIBRARY_PATH"] = os.pathsep.join([
      BUILD_DIR.abspath,
      env["ENV"].get("DYLD_FALLBACK_LIBRARY_PATH", "")])
else:
  env["ENV"]["LD_LIBRARY_PATH"] = os.pathsep.join([
      BUILD_DIR.abspath, env["ENV"].get("LD_LIBRARY_PATH", "")])

# Add the pysgpp package path to the environment
#########################################################################

if env["PLATFORM"] == "win32":
  # try to import pysgpp to detect an already existing installation, which
  # could cause trouble
  try:
    import pysgpp
    Helper.printWarning("An existing installation of pysgpp was detected. "
                        "To get rid of this warning remove the pysgpp package "
                        "from your local Python installation.")
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
  env["ENV"]["PYTHONPATH"] = os.pathsep.join([
      pysgppTempFolder,
                                              env["ENV"].get("PYTHONPATH", "")])
else:
  env["ENV"]["PYTHONPATH"] = os.pathsep.join([
      PYSGPP_PACKAGE_PATH.abspath,
      PYSGPP_BUILD_PATH.abspath,
      env["ENV"].get("PYTHONPATH", "")])

# Style checker
#########################################################################

def checkStyleAction(target, source, env):
  def run(args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # wait for termination and get output on stdout and stderr
    stdout, stderr = p.communicate()
    # in Python 3.x, communicate returns bytes
    if sys.version_info >= (3, 0):
      stdout, stderr = stdout.decode(), stderr.decode()
    return stdout, stderr

  sourcePath = source[0].abspath
  stdouts, stderrs = zip(*[run(args) for args in [
      ["python", "tools/cpplint.py", "--extensions=cpp,hpp",
       "--linelength=100", sourcePath],
      ["python", "tools/check_copyright_banners.py", sourcePath],
      ["python", "tools/check_includes.py", sourcePath]]])
  stdout, stderr = "\n".join(stdouts), "\n".join(stderrs)

  # cpplint prints on stderr
  for line in stderr.splitlines():
    # skip status lines, empty lines, and some warning types
    if ("Done processing " in line) or \
        ("Total errors found: " in line) or \
        ("Is this a non-const reference? " +
         "If so, make const or use a pointer:" in line) or \
        ("Consider using rand_r(...) instead of rand(...) for " +
         "improved thread safety." in line) or \
        (line == ""):
      pass
    else:
      parts = line.split(":  ")
      location = parts[0]
      message = ":  ".join(parts[1:])
      if message != "":
        print(location + ": warning: " + message)
      else:
        # occurs when cpplint excludes file via CPPLINT.cfg
        print(location)
  # touch file without writing anything
  # (to indicate for the next run of SCons that we already checked this file)
  with open(target[0].abspath, "w"): pass

env.Export("checkStyleAction")

# Custom builders for running Python/Boost tests and examples
#########################################################################

if env["RUN_PYTHON_TESTS"]:
  if env["SG_PYTHON"]:
    builder = Builder(action="python3 $SOURCE", chdir=0)
    env.Append(BUILDERS={"Test" : builder})
    builder = Builder(action="python3 $SOURCE")
    env.Append(BUILDERS={"SimpleTest" : builder})
  else:
    Helper.printWarning("Python tests disabled because SG_PYTHON is disabled.")

if env["COMPILE_BOOST_TESTS"]:
  builder = Builder(action="./$SOURCE --log_level=test_suite")
  env.Append(BUILDERS={"BoostTest" : builder})

if env["RUN_CPP_EXAMPLES"]:
  builder = Builder(action="./${SOURCE.file}", chdir=1)
  env.Append(BUILDERS={"CppExample" : builder})

if env["RUN_PYTHON_EXAMPLES"]:
  if env["SG_PYTHON"]:
    builder = Builder(action="python3 ${SOURCE.file}", chdir=1)
    env.Append(BUILDERS={"PythonExample" : builder})
  else:
    Helper.printWarning("Running of Python examples disabled "
                        "because SG_PYTHON is disabled.")

# Building the modules
#########################################################################

libraryTargetList = []
pythonTestTargetList = []
boostTestTargetList = []
boostTestRunTargetList = []
exampleTargetList = []
cppTestRunTargetList = []
pythonTestRunTargetList = []
pydocTargetList = []
headerSourceList = []
headerDestList = []
env.Export("libraryTargetList")
env.Export("pythonTestTargetList")
env.Export("boostTestTargetList")
env.Export("boostTestRunTargetList")
env.Export("exampleTargetList")
env.Export("cppTestRunTargetList")
env.Export("pythonTestRunTargetList")
env.Export("pydocTargetList")
env.Export("headerSourceList")
env.Export("headerDestList")

# compile selected modules
flattenedDependencyGraph = []

for moduleFolder in moduleFolders:
  if (not env["SG_" + moduleFolder.upper()]) or env.GetOption("help"):
    continue

  Helper.printInfo("Preparing to build module: {}".format(moduleFolder))
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

# compile matsgpp
if env["SG_MATLAB"]:
  env.SConscript("#/matsgpp/SConscript", {"env": env, "moduleName": "matsgpp"})

if ((not env["SG_PYTHON"]) and (not env["SG_JAVA"]) and (not env["SG_MATLAB"]) and
    (len(flattenedDependencyGraph) == 0)):
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

if env["RUN_CPP_EXAMPLES"]:
  env.Depends(cppTestRunTargetList, finalStepDependencies)
  finalStepDependencies.append(cppTestRunTargetList)
  env.SideEffect("sideEffectFinalSteps", cppTestRunTargetList)

if env["RUN_PYTHON_EXAMPLES"]:
  env.Depends(pythonTestRunTargetList, finalStepDependencies)
  finalStepDependencies.append(pythonTestRunTargetList)
  env.SideEffect("sideEffectFinalSteps", pythonTestRunTargetList)

# System-wide installation
#########################################################################

installLibSGpp = env.Alias("install-lib-sgpp",
                           env.Install(os.path.join(env.get("LIBDIR")),
                                       libraryTargetList))

headerFinalDestList = []
for headerDest in headerDestList:
  headerFinalDestList.append(os.path.join(env.get("INCLUDEDIR"), os.path.relpath(headerDest).split(os.sep, 2)[2]))

installIncSGpp = env.Alias("install-inc-sgpp",
                           env.InstallAs(headerFinalDestList, headerSourceList))

env.Alias("install", [installLibSGpp, installIncSGpp])

# Doxygen
#########################################################################

doxygen = env.Command("doc/xml/index.xml", "Doxyfile", "doxygen $SOURCE")
env.AddPostAction(doxygen, DoxygenHelper.patchNavtree)
env.Alias("doxygen", doxygen)
# SCons doesn't know the *.doxy dependencies of the doxygen target
# ==> always consider out-of-date with AlwaysBuild
env.AlwaysBuild(doxygen)

# Things to be cleaned
#########################################################################

# build + log files
env.Clean("clean", ["lib/", "jsgpp/java/", "config.log", "build.log"])
# Doxygen stuff
env.Clean("clean", ["Doxyfile", "doxygen_warnings.log", "doc/html/", "doc/xml/",
                    "base/doc/doxygen/modules.doxy"])

for module in moduleFolders:
  # PYDOC stuff
  env.Clean("clean", [module + "/Doxyfile", module + "/doc/doxygen_warnings.log",
                      module + "/doc/xml/"])

# Default targets
#########################################################################

finalMessagePrinter.sgppBuildPath = BUILD_DIR.abspath

finalMessagePrinter.pysgppPackagePath = (
    PYSGPP_PACKAGE_PATH.abspath + ":" + PYSGPP_BUILD_PATH.abspath)

if env.GetOption("help"):
  finalMessagePrinter.disable()
elif env.GetOption("clean"):
  env.Default(finalStepDependencies + ["clean"])
  finalMessagePrinter.disable()
else:
  env.Default(finalStepDependencies)
  if "doxygen" in BUILD_TARGETS:
    finalMessagePrinter.disable()
  elif not env["PRINT_INSTRUCTIONS"]:
    finalMessagePrinter.disable()

# dirty fix for Debian bug #893740
# (https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=893740),
# they seem to have reintroduced the Python-2-only syntax "dict.has_key"
# (instead of "in") in SCons/Script/Main.py, line 1111,
# occurs when cleaning, i.e., `scons -c`, while using Python 3.x for SCons
if not hasattr(os.environ, "has_key"):
  os.environ.has_key = (lambda x: x in os.environ)

# Save build variables to file
@atexit.register
def say_bye():
  vars.Save("buildVars.out", env)
