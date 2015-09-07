# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import glob
import SCons
import os
import SGppConfigure
from SCons.Script.SConscript import SConsEnvironment
import warnings

from Helper import *
from posix import chdir

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 1)
print "Using SCons", SCons.__version__

sconsenv = SConsEnvironment()
scons_ver = sconsenv._get_major_minor_revision(SCons.__version__)
if scons_ver < (2, 3, 0):
  warnings.warn("You are using an older version of scons than we do!\nSGpp officially supports scons >= 2.3.0.\nThere are reports that it also compiles with scons >= 2.1.0.")

# to ignore folders containing a SConscript file, do the following:
# ignoreFolders = ["jsgpp"]
ignoreFolders = []

# find all modules
moduleFolders, languageSupport = getModules(ignoreFolders)

prepareDoxyfile(moduleFolders)

moduleNames = []
languageSupportNames = []

for moduleFolder in moduleFolders:
  moduleName = 'SG_' + moduleFolder.upper()
  moduleNames.append(moduleName)

for wrapper in languageSupport:
    if wrapper == "pysgpp":
      languageSupportNames.append('SG_PYTHON')
    elif wrapper == "jsgpp":
      languageSupportNames.append('SG_JAVA')

print "Available modules:", ", ".join(moduleNames)
print "Available language support:", ", ".join(languageSupportNames)

vars = Variables("custom.py")

# define the flags
vars.Add('CPPFLAGS', 'Set additional Flags, they are compiler-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
vars.Add('LINKFLAGS', 'Set additional Linker-flags, they are linker-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
# define the target
vars.Add('MARCH', 'Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', None)
vars.Add('TARGETCPU', "Sets the processor you are compiling for. 'default' means using gcc with standard configuration. Also available are: 'ICC', here Intel Compiler in version 11 or higher must be used", 'default')
vars.Add(BoolVariable('OPT', "Sets optimization on and off", False))
# for compiling on LRZ without errors: omit unit tests
vars.Add(BoolVariable('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False))
vars.Add(BoolVariable('SG_PYTHON', 'Build with python Support', 'SG_PYTHON' in languageSupportNames))
vars.Add(BoolVariable('SG_JAVA', 'Build with java Support', 'SG_JAVA' in languageSupportNames))

for moduleName in moduleNames:
  vars.Add(BoolVariable(moduleName, 'Build the module ' + moduleName, True))

vars.Add(BoolVariable('SSE3_FALLBACK', 'Tries to build as much as possible with SSE3 instead of AVX (intrinsics based functions won\'t work)', False))
vars.Add('OUTPUT_PATH', 'Path where built libraries are installed. Needs a trailing slash!', '')
vars.Add(BoolVariable('VERBOSE', 'Set output verbosity', False))
vars.Add('CMD_LOGFILE', 'Specifies a file to capture the build log', 'build.log')
vars.Add(BoolVariable('USE_OCL', 'Sets OpenCL enabled state (Only actually enabled if also the OpenCL environment variables are set)', False))
vars.Add(BoolVariable('COMPILE_BOOST_TESTS', 'Compile the test cases written using Boost Test.', True))
vars.Add(BoolVariable('RUN_BOOST_TESTS', 'Run the test cases written using Boost Test (only if COMPILE_BOOST_TESTS is true).', True))
vars.Add(BoolVariable('USE_DOUBLE_PRECISION', 'If disabled, SG++ will compile using single precision (floats).', True))

vars.Add(BoolVariable('USE_ARMADILLO', 'Sets if Armadillo should be used (only relevant for SGPP::optimization).', False))
vars.Add(BoolVariable('USE_EIGEN', 'Sets if Eigen should be used (only relevant for SGPP::optimization).', False))
vars.Add(BoolVariable('USE_GMMPP', 'Sets if Gmm++ should be used (only relevant for SGPP::optimization).', False))
vars.Add(BoolVariable('USE_UMFPACK', 'Sets if UMFPACK should be used (only relevant for SGPP::optimization).', False))

# initialize environment
env = Environment(variables=vars, ENV=os.environ)
if 'CXX' in ARGUMENTS:
  print "CXX: ", ARGUMENTS['CXX']
  env['CXX'] = ARGUMENTS['CXX']
if 'CC' in ARGUMENTS:
  print "CC: ", ARGUMENTS['CC']
  env['CC'] = ARGUMENTS['CC']
if 'CPPFLAGS' in ARGUMENTS:
  env['CPPFLAGS'] = ARGUMENTS['CPPFLAGS'].split(",")
if 'CFLAGS' in ARGUMENTS:
  env['CFLAGS'] = ARGUMENTS['CFLAGS']
if 'CPPDEFINES' in ARGUMENTS:
  defineDict = {}
  for define in ARGUMENTS['CPPDEFINES'].split(","):
    key, value = define.split("=")
    defineDict[key] = value
  env.AppendUnique(CPPDEFINES = defineDict)
  print env['CPPDEFINES']
env.Export('moduleNames')
env.Export('moduleFolders')

# Help Text
Help("""---------------------------------------------------------------------

Type: 'scons [parameters]' to build the libraries

There are compiler optimizations for different platforms which can be
specified via parameters.

Parameters can be set either by setting the corresponding environment
variables, or directly via the commandline, e.g.,
> scons VERBOSE=True
to enable verbose compilation.


Specifying the target, the following options are available:
    - default: using the gcc toolchain with OpenMP 2
    - ICC: using the ICC toolchain with OpenMP 3 with standard x86_64 options

For LRZ, please execute:
module load python
module load gcc/4.5

FOR LRZ and when using intel compiler, execute:
export LIBPATH=$LD_LIBRARY_PATH

---------------------------------------------------------------------

Parameters are:
""" +
vars.GenerateHelpText(env))

# adds trailing slashes were required and if not present
BUILD_DIR = Dir(os.path.join(env['OUTPUT_PATH'], 'lib', 'sgpp'))
Export('BUILD_DIR')
PYSGPP_PACKAGE_PATH = Dir(os.path.join(env['OUTPUT_PATH'], 'lib'))
Export('PYSGPP_PACKAGE_PATH')
PYSGPP_BUILD_PATH = Dir(os.path.join(PYSGPP_PACKAGE_PATH.abspath, 'pysgpp'))
Export('PYSGPP_BUILD_PATH')
JSGPP_BUILD_PATH = Dir(os.path.join(env['OUTPUT_PATH'], 'lib', 'jsgpp'))
Export('JSGPP_BUILD_PATH')
EXAMPLE_DIR = Dir(os.path.join(env['OUTPUT_PATH'], 'bin', 'examples'))
Export('EXAMPLE_DIR')

# no checks if clean:
if not env.GetOption('clean'):
    SGppConfigure.doConfigure(env, moduleFolders, languageSupport)

# add #/lib/sgpp to LIBPATH
# (to add corresponding -L... flags to linker calls)
env.Append(LIBPATH=[BUILD_DIR])

# add C++ defines for all modules
cppdefines = []
for module in moduleNames:
    if env[module]:
        cppdefines.append(module)
env.Append(CPPDEFINES=cppdefines)

# environement setup finished, export environment
Export('env')

env.Append(CPPPATH=['#/tools'])
config = env.Configure()
Export('config')
# set up paths (Only Tested on Ubuntu!)
env["ENV"]["LD_LIBRARY_PATH"] = ":".join([
    env["ENV"].get("LD_LIBRARY_PATH", ""),
    BUILD_DIR.abspath])
env["ENV"]["PYTHONPATH"] = ":".join([
    env["ENV"].get("PYTHONPATH", ""),
    PYSGPP_PACKAGE_PATH.abspath])

# add custom builder to trigger the unittests after the build and to enable a special import test
if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
    # do the actual thing
    builder = Builder(action="python $SOURCE.file", chdir=1)
    env.Append(BUILDERS={'Test' : builder})
    builder = Builder(action="python $SOURCE")
    env.Append(BUILDERS={'SimpleTest' : builder})

if env['COMPILE_BOOST_TESTS']:
    builder = Builder(action="./$SOURCE")
    env.Append(BUILDERS={'BoostTest' : builder})



libraryTargetList = []
installTargetList = []
testTargetList = []
boostTestTargetList = []
exampleTargetList = []
env.Export('libraryTargetList')
env.Export('installTargetList')
env.Export('testTargetList')
env.Export('boostTestTargetList')
env.Export('exampleTargetList')

# compile selected modules
for moduleFolder in moduleFolders:
  if not env['SG_' + moduleFolder.upper()]:
    continue
  print "Preparing to build module: ", moduleFolder
  # SConscript('src/sgpp/SConscript' + moduleFolder, variant_dir='#/tmp/build/', duplicate=0)
  env.SConscript('#/' + moduleFolder + '/SConscript', {'env': env, 'moduleName': moduleFolder})

if env['SG_PYTHON']:
  env.SConscript('#/pysgpp/SConscript', {'env': env, 'moduleName': "pysgpp"})

if env['SG_JAVA']:
  env.SConscript('#/jsgpp/SConscript', {'env': env, 'moduleName': "jsgpp"})

# Python tests
#########################################################################

# execute first test after installing the last module
dependencies = [installTargetList]
separator = 70 * "-"

def printRunningPythonTests(target, source, env):
  print "\n" + separator + "\nRunning Python tests...\n" + separator

if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
  dependencies.append(env.Command('printRunningPythonTests', [], printRunningPythonTests))

  # serialize tests and move them at the end of the build
  for testTarget in testTargetList:
    env.Requires(testTarget, installTargetList)
    dependencies.append(testTarget)

# Boost tests
#########################################################################

def printRunningBoostTests(target, source, env):
  print "\n" + separator + "\nRunning Boost tests...\n" + separator

if env['COMPILE_BOOST_TESTS'] and env['RUN_BOOST_TESTS']:
  dependencies.append(env.Command('printRunningBoostTests', [], printRunningBoostTests))

  for testTarget in boostTestTargetList:
    env.Requires(testTarget, installTargetList)
    dependencies.append(testTarget)

# Examples
#########################################################################

def printLinkingExamples(target, source, env):
  print "\n" + separator + "\nLinking examples...\n" + separator

dependencies.append(env.Command('printLinkingExamples', [], printLinkingExamples))

for exampleTarget in exampleTargetList:
  env.Requires(exampleTarget, installTargetList)
  dependencies.append(exampleTarget)

# Final output
#########################################################################

def printFinished(target, source, env):
  import string
  with open("INSTRUCTIONS") as f:
    instructionsTemplate = string.Template(f.read())
  print
  print instructionsTemplate.safe_substitute(SGPP_BUILD_PATH=BUILD_DIR.abspath,
                                             PYSGPP_PACKAGE_PATH=PYSGPP_PACKAGE_PATH.abspath)

dependencies.append(env.Command('printFinished', [], printFinished))

# necessary to enforce an order on the final steps of the building of the wrapper
for i in range(len(dependencies) - 1):
  env.Depends(dependencies[i + 1], dependencies[i])
