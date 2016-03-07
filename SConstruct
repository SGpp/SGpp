# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import glob
import os
import subprocess
import SCons
from SCons.Script.SConscript import SConsEnvironment

import Helper
import SGppConfigure

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 1)
print "Using SCons", SCons.__version__

scons_ver = SConsEnvironment()._get_major_minor_revision(SCons.__version__)
if scons_ver < (2, 3, 0):
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
vars.Add('CPPFLAGS', 'Set additional Flags, they are compiler-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=Helper.multiParamConverter)
vars.Add('LINKFLAGS', 'Set additional Linker-flags, they are linker-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=Helper.multiParamConverter)
vars.Add('CPPPATH', 'Set paths where to look for additional headers', '')
vars.Add('LIBPATH', 'Set paths where to look for additional libs', '')
# define the target
vars.Add('ARCH', 'Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', 'sse3')
vars.Add('COMPILER', "Sets the processor you are compiling for. 'gnu' means using gcc with standard configuration. Also available are: 'intel', here Intel Compiler in version 11 or higher must be used", 'gnu')
vars.Add(BoolVariable('OPT', "Sets optimization on and off", False))
# for compiling on LRZ without errors: omit unit tests
vars.Add(BoolVariable('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False))
vars.Add(BoolVariable('SG_PYTHON', 'Build with python Support', 'SG_PYTHON' in languageSupportNames))
vars.Add(BoolVariable('PYDOC', 'Build python wrapper with comments', 'SG_PYTHON' in languageSupportNames))
vars.Add(BoolVariable('SG_JAVA', 'Build with java Support', 'SG_JAVA' in languageSupportNames))


for moduleName in moduleNames:
  vars.Add(BoolVariable(moduleName, 'Build the module ' + moduleName, True))

vars.Add('PREFIX', 'Path where built libraries, headers and other built files are installed (architecture independant). Needs a trailing slash!', '/usr/local')
vars.Add('EPREFIX', 'Path of architecture dependant files.' )
vars.Add('LIBDIR', 'Path of libraries.' )
vars.Add('INCLUDEDIR', 'Path of header files.' )
vars.Add(BoolVariable('VERBOSE', 'Set output verbosity', False))
vars.Add('CMD_LOGFILE', 'Specifies a file to capture the build log', 'build.log')
vars.Add(BoolVariable('USE_OCL', 'Sets OpenCL enabled state (Only actually enabled if also the OpenCL environment variables are set)', False))
vars.Add('OCL_INCLUDE_PATH', 'Specifies the location of the OpenCL header files (parent directory of "CL/").')
vars.Add('OCL_LIBRARY_PATH', 'Specifies the location of the OpenCL library.')
vars.Add('BOOST_INCLUDE_PATH', 'Specifies the location of the boost header files.', '/usr/include')
vars.Add('BOOST_LIBRARY_PATH', 'Specifies the location of the boost library.', '/usr/lib/x86_64-linux-gnu')
vars.Add(BoolVariable('COMPILE_BOOST_TESTS', 'Compile the test cases written using Boost Test.', True))
vars.Add(BoolVariable('COMPILE_BOOST_PERFORMANCE_TESTS', 'Compile the performance tests written using Boost Test. Currently only buildable with OpenCL enabled', False))
vars.Add(BoolVariable('RUN_BOOST_TESTS', 'Run the test cases written using Boost Test (only if COMPILE_BOOST_TESTS is true).', True))
vars.Add(BoolVariable('RUN_CPPLINT', 'Check compliance to Google\'s style guide using cpplint.', True))

vars.Add(BoolVariable('USE_ARMADILLO', 'Sets if Armadillo should be used (only relevant for sgpp::optimization).', False))
vars.Add(BoolVariable('USE_EIGEN', 'Sets if Eigen should be used (only relevant for sgpp::optimization).', False))
vars.Add(BoolVariable('USE_GMMPP', 'Sets if Gmm++ should be used (only relevant for sgpp::optimization).', False))
vars.Add(BoolVariable('USE_UMFPACK', 'Sets if UMFPACK should be used (only relevant for sgpp::optimization).', False))
vars.Add(BoolVariable('USE_STATICLIB', 'Sets if a static library should be built.', False))
vars.Add(BoolVariable('PRINT_INSTRUCTIONS', 'Print instruction for installing SG++.', True))

# create temporary environment to check which system and compiler we should use
# (the Environment call without "tools=[]" crashes with MinGW,
# so we do it like that)
env = Environment(variables=vars, ENV=os.environ, tools=[])

if (env['PLATFORM'].lower() == 'win32') and \
   (env['COMPILER'].lower() == 'gnu'):
  # MinGW: use gcc toolschain
  tools = ['gnulink', 'gcc', 'g++', 'gas', 'ar', 'swig']
else:
  # otherwise: use default toolchain
  tools = ['default']

# initialize environment
env = Environment(variables=vars, ENV=os.environ, tools=tools)

env['EPREFIX'] = env.get( 'EPREFIX', env['PREFIX'] )
env['LIBDIR'] = env.get( 'LIBDIR', os.path.join( env['EPREFIX'], "lib" ) )
env['INCLUDEDIR'] = env.get( 'INCLUDEDIR', os.path.join( env['PREFIX'], "include" ) )

# no docu if clean:
if not env.GetOption('clean'):
  Helper.prepareDoxyfile(moduleFolders)

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
env.AppendUnique(CPPDEFINES = {})
if 'CPPDEFINES' in ARGUMENTS:
  for define in ARGUMENTS['CPPDEFINES'].split(","):
    key, value = define.split("=")
    env['CPPDEFINES'][key] = value

if 'CPPPATH' in ARGUMENTS:
    env['CPPPATH'] = ARGUMENTS['CPPPATH'].split(",")
if 'LIBPATH' in ARGUMENTS:
    env['LIBPATH'] = ARGUMENTS['LIBPATH'].split(",")

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
BUILD_DIR = Dir(os.path.join('lib', 'sgpp'))
Export('BUILD_DIR')
PYSGPP_PACKAGE_PATH = Dir(os.path.join('lib'))
Export('PYSGPP_PACKAGE_PATH')
PYSGPP_BUILD_PATH = Dir(os.path.join(PYSGPP_PACKAGE_PATH.abspath, 'pysgpp'))
Export('PYSGPP_BUILD_PATH')
JSGPP_BUILD_PATH = Dir(os.path.join('lib', 'jsgpp'))
Export('JSGPP_BUILD_PATH')
EXAMPLE_DIR = Dir(os.path.join('bin', 'examples'))
Export('EXAMPLE_DIR')

# no checks if clean:
if not env.GetOption('clean'):
  SGppConfigure.doConfigure(env, moduleFolders, languageSupport)

# fix for "command line too long" errors on MinGW
# (from https://bitbucket.org/scons/scons/wiki/LongCmdLinesOnWin32)
if env['PLATFORM'] == 'win32':
  Helper.set_win32_spawn(env)

# add #/lib/sgpp to LIBPATH
# (to add corresponding -L... flags to linker calls)
env.Append(LIBPATH=[BUILD_DIR])

# # add C++ defines for all modules
for module in moduleNames:
    if env[module]:
        env['CPPDEFINES'][module] = '1'

# environement setup finished, export environment
Export('env')

env.Append(CPPPATH=['#/tools'])
config = env.Configure()
Export('config')

# update PATH under win32/LD_LIBRARY_PATH otherwise
# to add BUILD_DIR (so we can run the Boost tests)
if env['PLATFORM'] == 'win32':
    env["ENV"]["PATH"] = os.pathsep.join([
        env["ENV"].get("PATH", ""),
        BUILD_DIR.abspath])

    # also add the Boost library path to the PATH
    # so that the Boost test *.dll can be found when running the tests
    if env['RUN_BOOST_TESTS']:
      env["ENV"]["PATH"] = os.pathsep.join([
          env["ENV"].get("PATH", ""),
          env["BOOST_LIBRARY_PATH"]])

# Mac OS X doens't use LD_LIBRARY_PATH 
elif env['PLATFORM'] == 'darwin':
    env["ENV"]["DYLD_FALLBACK_LIBRARY_PATH"] = os.pathsep.join([
        env["ENV"].get("DYLD_FALLBACK_LIBRARY_PATH", ""),
        BUILD_DIR.abspath])

else:
    env["ENV"]["LD_LIBRARY_PATH"] = os.pathsep.join([
        env["ENV"].get("LD_LIBRARY_PATH", ""),
        BUILD_DIR.abspath])

# -------------------------------------------------------------------------
# add the pysgpp package path to the environment
if env['PLATFORM'] == 'win32':
  # try to import pysgpp to detect an already existing installation, which
  # could cause trouble
  try:
    import pysgpp
    print "Warning: more than one installations of pysgpp are detected. To get rid of this warning remove the pysgpp package from your local python installation."
  except:
    pass

  # get a temporary folders
  import tempfile, uuid
  # get temp directory
  pysgppTempFolder = os.path.join(tempfile.gettempdir(),
                           "site-pyspp-%s" % str(uuid.uuid4()))
  # create temp folder
  os.makedirs(pysgppTempFolder)

  # add it to the build python path
  env["ENV"]["PYTHONPATH"] = os.pathsep.join([pysgppTempFolder,
                                              env["ENV"].get("PYTHONPATH", "")])
else:
  env["ENV"]["PYTHONPATH"] = os.pathsep.join([env["ENV"].get("PYTHONPATH", ""),
                                              PYSGPP_PACKAGE_PATH.abspath])
# -------------------------------------------------------------------------

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

env.Export('lintAction')

# -------------------------------------------------------------------------

# add custom builder to trigger the unittests after the build and to enable a special import test
if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
    # do the actual thing
    builder = Builder(action="python $SOURCE", chdir=0)
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
pydocTargetList = []
headerSourceList = []
headerDestList = []
env.Export('libraryTargetList')
env.Export('installTargetList')
env.Export('testTargetList')
env.Export('boostTestTargetList')
env.Export('exampleTargetList')
env.Export('pydocTargetList')
env.Export('headerSourceList')
env.Export('headerDestList')

# compile selected modules
flattenedDependencyGraph = []
for moduleFolder in moduleFolders:
  if not env['SG_' + moduleFolder.upper()]:
    continue
  print "Preparing to build module: ", moduleFolder
  # SConscript('src/sgpp/SConscript' + moduleFolder, variant_dir='#/tmp/build/', duplicate=0)
  env.SConscript('#/' + moduleFolder + '/SConscript', {'env': env, 'moduleName': moduleFolder})

  # add the dependencies of the current module to the overall dependency graph
  Import("moduleDependencies")
  Import("libname")
  flattenedDependencyGraph = Helper.flatDependencyGraph(
      [libname] + moduleDependencies, flattenedDependencyGraph)

Export('flattenedDependencyGraph')

if env['PYDOC'] and env['SG_PYTHON']:
  data = open('moduleDoxy', 'r').read()
  for module in moduleFolders:
    if not env['SG_' + module.upper()]:
      continue

    if env.GetOption("clean"):
      if os.path.exists(os.path.join(module, 'Doxyfile')):
        os.remove(os.path.join(module, 'Doxyfile'))
      doxypath = os.path.join(module, 'doc/xml/')
      if os.path.exists(doxypath):
        for file in os.listdir(doxypath):
          os.remove(os.path.join(doxypath, file))
    else:
      with open(os.path.join(module, 'Doxyfile'), 'w') as doxyFile:
        doxyFile.write(data.replace('$modname', module).replace('$quiet', 'YES'))

    doxy_env = env.Clone()

    doxygen = doxy_env.Command(os.path.join(module, 'doc/xml/index.xml'), '', 'doxygen ' + os.path.join(module, 'Doxyfile'))

    doxy2swig_command = "python pysgpp/doxy2swig.py -o -c -q $SOURCE $TARGET"
    doxy2swig = doxy_env.Command(os.path.join('pysgpp', module + '_doc.i'), doxygen, doxy2swig_command)

    for root, dirs, files in os.walk(os.path.join(module, 'src')):
      for file in files:
        if 'cpp' in file or 'hpp' in file:
          doxy_env.Depends(doxygen, os.path.join(root, file))
          doxy_env.Depends(doxy2swig, os.path.join(root, file))
    pydocTargetList.append(doxy2swig)

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

def installPythonLibToTmp(target, source, env):
  # prepare python package for unit testing
  import sys, os, subprocess

  # get temp directory
  pysgppTempFolder = source[0].get_string(0)

  # install python interface to tmp directory
  p = subprocess.call(["python", "setup.py",
                       "--quiet",
                       "install", "--install-lib=%s" % pysgppTempFolder])
  if p != 0:
      print "Error: installing python package to the temporary folder '%s' failed; I can not run the python unit tests automatically." % pysgppTempFolder
      exit(-1)

if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
  if env['PLATFORM'] == 'win32':
    # install the python library to that temporary folder
    dependencies.append(env.Command('installPythonLibToTmp', [pysgppTempFolder], installPythonLibToTmp))
    

  # print message that python tests are about to start
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
  if env['PLATFORM'] in ['cygwin', 'win32']:
    filename = "INSTRUCTIONS_WINDOWS"
  elif env['PLATFORM'] == 'darwin' :
    filename = "INSTRUCTIONS_MAC"
  else:
    filename = "INSTRUCTIONS"
 
  with open(filename) as f:
    instructionsTemplate = string.Template(f.read())
    print
    print instructionsTemplate.safe_substitute(SGPP_BUILD_PATH=BUILD_DIR.abspath,
                                               PYSGPP_PACKAGE_PATH=PYSGPP_PACKAGE_PATH.abspath)
 
if env["PRINT_INSTRUCTIONS"]:
    dependencies.append(env.Command('printFinished', [], printFinished))

# necessary to enforce an order on the final steps of the building of the wrapper
for i in range(len(dependencies) - 1):
  env.Depends(dependencies[i + 1], dependencies[i])
 
# Stuff needed for system install
env.Clean("distclean",
  [
    "config.log",
  ])
# TODO(killian): the next line messes up the dependency tracking, seems to not be required, please check (reported by David)
#Default(libraryTargetList, dependencies)
    
ils = env.Alias('install-lib-sgpp', Install(os.path.join( env.get('LIBDIR'), 'sgpp'), libraryTargetList))
   
headerFinalDestList = []
for headerDest in headerDestList:
  headerFinalDestList.append(os.path.join( env.get('INCLUDEDIR'), headerDest))
   
iis = env.Alias('install-inc-sgpp', InstallAs(headerFinalDestList, headerSourceList))
   
env.Alias('install', [ils, iis])
