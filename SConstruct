# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import glob
import SCons
import fnmatch
import os
import SGppConfigure

from Helper import *

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 0)
print "Using SCons", SCons.__version__

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

print moduleFolders
print moduleNames

print languageSupport
print languageSupportNames

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

# initialize environment
env = Environment(variables=vars, ENV=os.environ)
if 'CXX' in ARGUMENTS:
  print "CXX: ", ARGUMENTS['CXX']
  env['CXX'] = ARGUMENTS['CXX']
if 'CC' in ARGUMENTS:
  print "CC: ", ARGUMENTS['CC']
  env['CC'] = ARGUMENTS['CC']
if 'CPPFLAGS' in ARGUMENTS:
  env['CPPFLAGS'] = ARGUMENTS['CPPFLAGS'].split(" ")
if 'CFLAGS' in ARGUMENTS:
  env['CFLAGS'] = ARGUMENTS['CFLAGS']
if 'CPPDEFINES' in ARGUMENTS:
  defineDict = {}
  for define in ARGUMENTS['CPPDEFINES'].split(" "):
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
PYSGPP_BUILD_PATH = Dir(os.path.join(env['OUTPUT_PATH'], 'lib', 'pysgpp'))
Export('PYSGPP_BUILD_PATH')
JSGPP_BUILD_PATH = Dir(os.path.join(env['OUTPUT_PATH'], 'lib', 'jsgpp'))
Export('JSGPP_BUILD_PATH')
EXAMPLE_DIR = Dir(os.path.join(env['OUTPUT_PATH'], 'bin', 'examples'))
Export('EXAMPLE_DIR')

# no checks if clean:
if not env.GetOption('clean'):
    SGppConfigure.doConfigure(env, moduleFolders, languageSupport)

# add #/lib/sgpp and #/lib/alglib to LIBPATH
# (to add corresponding -L... flags to linker calls)
env.Append(LIBPATH=[BUILD_DIR,
                    Dir(os.path.join(env['OUTPUT_PATH'], 'lib', 'alglib'))])

# add C++ defines for all modules
cppdefines = []
for module in moduleNames:
    cppdefines.append(module)
env.Append(CPPDEFINES=cppdefines)

# environement setup finished, export environment
Export('env')

# Install alglib
libalglib, alglibstatic = SConscript(os.path.join('tools', 'SConscriptAlglib'),
                                      variant_dir=os.path.join('tmp', 'build_alglib'),
                                      duplicate=0)
alglibinst = env.Install(os.path.join(env['OUTPUT_PATH'], 'lib', 'alglib'),
                         [libalglib, alglibstatic])
env.Depends(os.path.join("#", BUILD_DIR.path, "libsgppbase.so"), alglibinst)

env.Append(CPPPATH=['#/tools'])

# add custom builder to trigger the unittests after the build and to enable a special import test
if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
    # run tests
    builder = Builder(action="python $SOURCE.file", chdir=1)
    env.Append(BUILDERS={'Test' : builder})
    builder = Builder(action="python $SOURCE")
    env.Append(BUILDERS={'SimpleTest' : builder})

libraryTargetList = []
installTargetList = []
testTargetList = []
exampleTargetList = []
env.Export('libraryTargetList')
env.Export('installTargetList')
env.Export('testTargetList')
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

#TODO: ask julian
if env['SG_JAVA']:
  env.SConscript('#/jsgpp/SConscript', {'env': env, 'moduleName': "jsgpp"})

# Unit tests
#########################################################################

# necessary to enforce an order on the final steps of the building of the wrapper
dependency = None
if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
  # serialize tests and move them at the end of the build
  for testTarget in testTargetList:
    env.Requires(testTarget, installTargetList)

    if dependency is None:
      #print testTarget, 'depends on nothing'
      dependency = testTarget
    else:
      #print testTarget, 'depends on', dependency
      env.Depends(testTarget, dependency)
      dependency = testTarget

for exampleTarget in exampleTargetList:
  env.Requires(exampleTarget, installTargetList)
  if dependency is None:
    dependency = exampleTarget
  else:
    env.Depends(exampleTarget, dependency)
    dependency = exampleTarget
