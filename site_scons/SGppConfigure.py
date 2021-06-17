# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


import distutils.sysconfig
import errno
import os
import re
import subprocess
import sys

import SCons.Script

import Helper

def getOutput(command):
  # redirect stderr to stdout
  try:
    output = subprocess.check_output(command, stderr=subprocess.STDOUT)
  except subprocess.CalledProcessError as e:
    output = e.output  
  # in Python 3.x, check_output returns bytes
  if sys.version_info >= (3, 0): output = output.decode()
  # strip trailing newlines
  output = output.rstrip("\r\n")
  return output

def doConfigure(env, moduleFolders, languageWrapperFolders):
  print("")
  print("Checking programs and libraries:")

  config = env.Configure(custom_tests={"CheckExec" : Helper.CheckExec,
                                       "CheckJNI" : Helper.CheckJNI,
                                       "CheckFlag" : Helper.CheckFlag,
                                       "CheckCompiler" : Helper.CheckCompiler,
                                       "CheckMKL" : Helper.CheckMklScalapack})

  # now set up all further environment settings that should never fail
  # compiler setup should be always after checking headers and flags,
  # as they can make the checks invalid, e.g., by setting "-Werror"

  if env["OPT"] == True:
    config.env.Append(CPPFLAGS=["-O3", "-g"])

  else:
    config.env.Append(CPPFLAGS=["-g", "-O0"])

  # make settings case-insensitive
  config.env["COMPILER"] = config.env["COMPILER"].lower()
  config.env["ARCH"] = config.env["ARCH"].lower()

  if config.env["COMPILER"] in ("gnu", "openmpi", "mpich"):
    configureGNUCompiler(config)
  elif config.env["COMPILER"] == "clang":
    configureClangCompiler(config)
  elif config.env["COMPILER"] in ("intel", "intel.mpi"):
    configureIntelCompiler(config)
  else:
    Helper.printErrorAndExit("You must specify a valid value for Compiler.",
                             "Available configurations are:",
                             "gnu, clang, intel, openmpi, mpich, intel.mpi")

  if config.env["COMPILER"] in ("openmpi", "mpich", "intel.mpi"):
    config.env["CPPDEFINES"]["USE_MPI"] = "1"
    config.env["USE_MPI"] = True # tells scons to build MPI related examples and operations
  else:
    config.env["USE_MPI"] = False

  # special treatment for different platforms
  if config.env["PLATFORM"] == "darwin":
    # the "-undefined dynamic_lookup"-switch is required to actually build a shared library
    # in OSX. "-dynamiclib" alone results in static linking of
    # all further dependent shared libraries
    # beware: if symbols are missing that are actually required
    # (because the symbols don't reside in a shared library),
    # there will be no error during compilation
    # also for the python binding, the library must be suffixed with '*.so' even
    # though it is a dynamiclib and not a bundle (see SConscript in src/pysgpp)
    config.env.AppendUnique(LINKFLAGS=["-undefined", "dynamic_lookup"])
    # The GNU assembler (GAS) is not supported in Mac OS X.
    # A solution that fixed this problem is by adding -Wa,-q to the compiler flags.
    # From the man pages for as (version 1.38):
    # -q Use the clang(1) integrated assembler instead of the GNU based system assembler.
    # Note that the CPPFLAG is exactly "-Wa,-q", where -Wa passes flags to the assembler and
    # -q is the relevant flag to make it use integrated assembler
    if config.env["COMPILER"] == "gcc":
      config.env.AppendUnique(CPPFLAGS=["-Wa,-q"])
    config.env.AppendUnique(CPPPATH="/usr/local/include")
    config.env.AppendUnique(LIBPATH="/usr/local/lib")
    config.env["SHLIBSUFFIX"] = ".dylib"
  elif config.env["PLATFORM"] == "cygwin":
    # required to find the static libraries compiled before the shared libraries
    # the static libraries are required as the linker on windows cannot ignore undefined symbols
    # (as is done on linux automatically and is done on OSX with the settings above)
    #env.Append(LIBPATH=[BUILD_DIR])
    pass

  # will lead to a warning on cygwin and win32
  # ("all code is position independent")
  if config.env["PLATFORM"] not in ["cygwin", "win32"]:
    config.env.AppendUnique(CPPFLAGS=["-fPIC"])

  # setup the include base folder
  # env.Append(CPPPATH=["#/src/sgpp"])
  for moduleFolder in moduleFolders:
    if (moduleFolder in languageWrapperFolders) or (not config.env["SG_" + moduleFolder.upper()]):
      continue
    config.env.AppendUnique(CPPPATH=["#/" + moduleFolder + "/src/"])

  # detour compiler output
  config.env["PRINT_CMD_LINE_FUNC"] = Helper.printCommand

  checkCpp11(config)
  if "doxygen" in SCons.Script.BUILD_TARGETS:
    checkDoxygen(config)
    checkDot(config)
  checkOpenCL(config)
  detectGSL(config)
  detectZlib(config)
  detectScaLAPACK(config)
  checkDAKOTA(config)
  checkCGAL(config)
  checkBoostTests(config)
  checkSWIG(config)
  checkPython(config)
  checkJava(config)

  if config.env["USE_CUDA"] == True:
    config.env['CUDA_TOOLKIT_PATH'] = ''
    config.env['CUDA_SDK_PATH'] = ''
    config.env.Tool('cuda')
    # clean up the flags to forward
    # flagsToForward = [flag for flag in config.env["CPPFLAGS"] if flag not in ['-Wmissing-format-attribute', '']]
    # flagsToForward = " -Xcompiler " + (" -Xcompiler ".join(flagsToForward))
    # ensure same flags for host code
    config.env['NVCCFLAGS'] = "-ccbin " + config.env["CXX"] + " -std=c++11 -Xcompiler -fpic,-Wall "# + flagsToForward
    # config.env.AppendUnique(LIBPATH=['/usr/local.nfs/sw/cuda/cuda-7.5/'])

  if config.env["USE_HPX"]:
    hpxLibs = ["dl", "rt", "boost_chrono", "boost_date_time", "boost_filesystem", "boost_program_options", "boost_regex" ,
               "boost_system", "boost_thread", "boost_context", "boost_random", "boost_atomic", "tcmalloc_minimal", "hwloc"]
    if config.env["OPT"]:
      hpxLibs += ["hpx", "hpx_init", "hpx_iostreams"]
      if "HPX_RELEASE_LIBRARY_PATH" in config.env:
        config.env.AppendUnique(LIBPATH=config.env["HPX_RELEASE_LIBRARY_PATH"])
    else:
      if "HPX_DEBUG_LIBRARY_PATH" in env:
        config.env.AppendUnique(LIBPATH=config.env["HPX_DEBUG_LIBRARY_PATH"])
      hpxLibs += ["hpxd", "hpx_initd", "hpx_iostreamsd"]

    for lib in hpxLibs:
      if not config.CheckLib(lib, language="c++", autoadd=1):
        Helper.printErrorAndExit("lib" + lib + " not found, but required for HPX")
    config.env.AppendUnique(LIBS=hpxLibs)
    config.env["CPPDEFINES"]["USE_HPX"] = "1"


    if 'HPX_SHARED_INCLUDE_PATH' in config.env:
      config.env.AppendUnique(CPPPATH=config.env['HPX_SHARED_INCLUDE_PATH'])
    if config.env["OPT"]:
      config.env.AppendUnique(CPPDEFINES=["HPX_APPLICATION_EXPORTS", "HPX_ENABLE_ASSERT_HANDLER"]);
      if 'HPX_RELEASE_INCLUDE_PATH' in config.env:
        config.env.AppendUnique(CPPPATH=config.env['HPX_RELEASE_INCLUDE_PATH'])
    else:
      config.env.AppendUnique(CPPDEFINES=["HPX_DEBUG", "HPX_APPLICATION_EXPORTS", "HPX_ENABLE_ASSERT_HANDLER"]);
      if 'HPX_DEBUG_INCLUDE_PATH' in config.env:
        config.env.AppendUnique(CPPPATH=config.env['HPX_DEBUG_INCLUDE_PATH'])

    if not config.CheckCXXHeader("hpx/hpx_init.hpp"):
      Helper.printErrorAndExit("hpx/hpx_init.hpp not found, but required for HPX")
    if not config.CheckCXXHeader("hpx/include/actions.hpp"):
      Helper.printErrorAndExit("hpx/include/actions.hpp not found, but required for HPX")

  env = config.Finish()

  print("Configuration done.")
  print("")

def checkCpp11(config):
  # check C++11 support
  if not config.env['USE_HPX']:
    if not config.CheckFlag("-std=c++11"):
      Helper.printErrorAndExit("The compiler doesn't seem to support the C++11 standard. Abort!")
      Exit(1)

    config.env.AppendUnique(CPPFLAGS="-std=c++11")
  else:
    if not config.CheckFlag("-std=c++14"):
      Helper.printErrorAndExit("HPX requires a compiler that supports the C++14 standard. Abort!")
      Exit(1)

    config.env.AppendUnique(CPPFLAGS="-std=c++14")

def checkDoxygen(config):
  # check whether Doxygen installed
  if not config.CheckExec("doxygen"):
    Helper.printWarning("Doxygen cannot be found.",
                        "You will not be able to generate the documentation.",
                        "Check the PATH environment variable!")
  else:
    Helper.printInfo("Using Doxygen " + re.findall(
       r"[0-9.]*[0-9]+", getOutput(["doxygen", "--version"]))[0] + ".")

def checkDot(config):
  # check whether dot installed
  if not config.CheckExec("dot"):
    Helper.printWarning("dot (Graphviz) cannot be found.",
                        "The documentation might lack diagrams.",
                        "Check the PATH environment variable!")
  else:
    Helper.printInfo("Using " + getOutput(["dot", "-V"]) + ".")

def checkOpenCL(config):

  # OpenCL also need boost to build
  config.env.AppendUnique(CPPPATH=[config.env["BOOST_INCLUDE_PATH"]])
  config.env.AppendUnique(LIBPATH=[config.env["BOOST_LIBRARY_PATH"]])

  if config.env["USE_OCL"]:
    if "OCL_INCLUDE_PATH" in config.env["ENV"]:
      config.env.AppendUnique(CPPPATH=[config.env["ENV"]["OCL_INCLUDE_PATH"]])
    elif "OCL_INCLUDE_PATH" in config.env:
      config.env.AppendUnique(CPPPATH=[config.env["OCL_INCLUDE_PATH"]])
    else:
      Helper.printInfo("Trying to find the OpenCL header without the OCL_INCLUDE_PATH variable.")

    if not config.CheckCXXHeader("CL/cl.h"):
      Helper.printErrorAndExit("CL/cl.h not found, but required for OpenCL")

    if "OCL_LIBRARY_PATH" in config.env["ENV"]:
      config.env.AppendUnique(LIBPATH=[config.env["ENV"]["OCL_LIBRARY_PATH"]])
    elif "OCL_LIBRARY_PATH" in config.env:
      config.env.AppendUnique(LIBPATH=[config.env["OCL_LIBRARY_PATH"]])
    else:
      Helper.printInfo("Trying to find the libOpenCL library without the OCL_LIBRARY_PATH variable.")

    if not config.CheckLib("OpenCL", language="c++", autoadd=1):
      Helper.printErrorAndExit("libOpenCL not found, but required for OpenCL")

    if not config.CheckLib("boost_program_options", language="c++", autoadd=0):
      Helper.printErrorAndExit("libboost-program-options not found, but required for OpenCL",
                               "On debian-like system the package libboost-program-options-dev",
                               "can be installed to solve this issue.")

    if not config.CheckCXXHeader("zlib.h"):
      Helper.printErrorAndExit("zlib.h not found, but required for OpenCL",
                               "On debian-like system the package zlib1g-dev",
                               "can be installed to solve this issue.")

    if not config.CheckLib("libz", language="c++", autoadd=0):
      Helper.printErrorAndExit("libz not found, but required for OpenCL",
                               "On debian-like system the package zlib1g",
                               "can be installed to solve this issue.")

    config.env["CPPDEFINES"]["USE_OCL"] = "1"

def checkDAKOTA(config):
    if config.env["USE_DAKOTA"]:
        if not config.CheckCXXHeader("pecos_global_defs.hpp"):
            Helper.printErrorAndExit("pecos_global_defs.hpp not found, but required for PECOS. Consider setting the flag 'CPPPATH'.")

def checkCGAL(config):
    if config.env["USE_CGAL"]:
        if not config.CheckCXXHeader("CGAL/basic.h"):
            Helper.printErrorAndExit("CGAL/basic.h not found, but required for CGAL. Consider setting the flag 'CPPPATH'.")

def checkBoostTests(config):
  # Check the availability of the boost unit test dependencies
  if config.env["COMPILE_BOOST_TESTS"]:
    config.env.AppendUnique(CPPPATH=[config.env["BOOST_INCLUDE_PATH"]])
    config.env.AppendUnique(LIBPATH=[config.env["BOOST_LIBRARY_PATH"]])

    if not config.CheckHeader(os.path.join("boost", "test", "unit_test.hpp"), language="c++"):
      config.env["COMPILE_BOOST_TESTS"] = False
      Helper.printWarning("****************************************************",
                          "No Boost Unit Test Headers found. Omitting Boost unit tests.",
                          "Please install the corresponding package, e.g., on Ubuntu",
                          "> sudo apt-get install libboost-test-dev",
                          "****************************************************")

    if not config.CheckLib("boost_unit_test_framework", autoadd=0, language="c++"):
      config.env["COMPILE_BOOST_TESTS"] = False
      Helper.printWarning("****************************************************",
                          "No Boost Unit Test library found. Omitting Boost unit tests.",
                          "Please install the corresponding package, e.g., on Ubuntu",
                          "> sudo apt-get install libboost-test-dev",
                          "****************************************************")

def checkSWIG(config):
  if config.env["SG_PYTHON"] or config.env["SG_JAVA"]:
    # check whether swig installed
    if not config.CheckExec("swig"):
      Helper.printErrorAndExit("SWIG cannot be found, but required for SG_PYTHON.",
                               "Check the PATH environment variable!")

    # make sure swig version is new enough
    swigVersion = re.findall(
        r"[0-9.]*[0-9]+", getOutput(["swig", "-version"]))[0]

    swigVersionTuple = config.env._get_major_minor_revision(swigVersion)
    if swigVersionTuple < (3, 0, 4):
      Helper.printErrorAndExit("SWIG version too old! At least 3.0.4 required.")

    Helper.printInfo("Using SWIG {}".format(swigVersion))

def checkPython(config):
  if config.env["SG_PYTHON"]:
    if not python3_is_installed():
        raise Exception("Python 3 is required for SGpp python support!")
      
    pythonpath = getOutput(["python3", "-c",
          "import distutils.sysconfig; "
          "print(distutils.sysconfig.get_python_inc())"])
    package = "python3-dev"

    config.env.AppendUnique(CPPPATH=[pythonpath])
    Helper.printInfo("pythonpath = " + pythonpath)

    if not config.CheckCXXHeader("Python.h"):
      Helper.printErrorAndExit("Python.h not found, but required for SG_PYTHON.",
                               "Check path to Python include files:",
                               pythonpath,
                               "Hint: You might have to install the package " + package + ".")

    if not config.CheckCXXHeader("pyconfig.h"):
      Helper.printErrorAndExit("pyconfig.h not found, but required for SG_PYTHON.",
                               "Check path to Python include files:",
                               pythonpath,
                               "Hint: You might have to install the package " + package + ".")

    numpy_path=getOutput(["python3", "-c", "import numpy, os;"
    "print(os.path.join(os.path.split(numpy.__file__)[0], \"core\", \"include\"))"]) 
    if numpy_path.startswith("Traceback"):
      Helper.printWarning("Warning: Numpy doesn't seem to be installed.")
      if config.env["RUN_PYTHON_TESTS"]:
        Helper.printWarning("Python unit tests were disabled because numpy is not available.")
        config.env["RUN_PYTHON_TESTS"] = False
      if config.env["RUN_PYTHON_EXAMPLES"]:
        Helper.printWarning("Python examples were disabled because numpy is not available.")
        config.env["RUN_PYTHON_EXAMPLES"] = False
    else:
      config.env.AppendUnique(CPPPATH=[numpy_path])
      if not config.CheckCXXHeader(["Python.h", "pyconfig.h", "numpy/arrayobject.h"]):
        Helper.printWarning("Cannot find NumPy header files in " + str(numpy_path) + ".")
        if config.env["RUN_PYTHON_TESTS"]:
          config.env["RUN_PYTHON_TESTS"] = False
          Helper.printWarning("Python unit tests were disabled due to missing numpy development headers.")

    if getOutput(["python3", "-c", "import scipy; "]).startswith('Traceback'):
      Helper.printWarning("Warning: Scipy doesn't seem to be installed.")
  else:
    Helper.printInfo("Python extension (SG_PYTHON) not enabled.")

def python3_is_installed():
  try:
    subprocess.check_output(["python3", "--version"])
    return True
  except subprocess.CalledProcessError:
    return False
  except OSError as e:
    # file not found
    if e.errno == errno.ENOENT:
      return False
    else:
      raise

def checkJava(config):
  if config.env["SG_JAVA"]:
    # check for $JAVA_HOME; prepend to search path
    if os.environ.get("JAVA_HOME"):
      config.env.PrependENVPath("PATH", os.path.join(os.environ.get("JAVA_HOME"), "bin"))

    # check whether javac installed
    if not config.CheckExec("javac"):
      Helper.printErrorAndExit("javac cannot be found, but required by SG_JAVA.",
                               "Check the PATH environment variable!")

    # check whether java installed
    if not config.CheckExec("java"):
      Helper.printErrorAndExit("java cannot be found, but required by SG_JAVA.",
                               "Check the PATH environment variable!")

    # check for JNI headers
    if os.environ.get("JNI_CPPINCLUDE"):
      config.env.AppendUnique(CPPPATH=[os.environ.get("JNI_CPPINCLUDE")])
    if not config.CheckCXXHeader("jni.h"):
      # not found; try to find
      if not config.CheckJNI():
        Helper.printErrorAndExit(" jni.h not found."
                                 "Please set JAVA_HOME environment variable"
                                 "with $JAVA_HOME/bin/javac, $JAVA_HOME/include/jni.h"
                                 "or directly $JNI_CPPINCLUDE with $JNI_CPPINCLUDE/jni.h.")
  else:
    Helper.printInfo("Java support (SG_JAVA) not enabled.")

def configureGNUCompiler(config):

  if config.env["RUN_ON_HAZELHEN"]:
    if "CC" not in config.env.arguments: config.env["CC"] = "CC"
    if "CXX" not in config.env.arguments: config.env["CXX"] = "CC"
    config.env.Append(CPPPATH = [os.environ['BOOST_ROOT'] + '/include'])
    config.env.Append(LIBPATH = [os.environ['BOOST_ROOT'] + '/lib'])
    config.env.Append(CPPFLAGS=["-dynamic"])
    config.env.Append(LINKFLAGS=["-dynamic"])
  if config.env["COMPILER"] == "openmpi":
    if "CC" not in config.env.arguments: config.env["CC"] = "mpicc"
    if "LINK" not in config.env.arguments: config.env["LINK"] = "mpicxx"
    if "CXX" not in config.env.arguments: config.env["CXX"] = "mpicxx"
    Helper.printInfo("Using openmpi.")
    # openmpi specific fix according to: https://github.com/open-mpi/ompi/issues/5157
    config.env["CPPDEFINES"]["OMPI_SKIP_MPICXX"] = "1"
  elif config.env["COMPILER"] == "mpich":
    if config.env["CC"]:
      config.env.Append(CFLAGS=["-cc=" + config.env["CC"]])
    if config.env["CXX"]:
      config.env.Append(CPPFLAGS=["-cxx=" + config.env["CXX"]])
      config.env.Append(LINKFLAGS=["-cxx=" + config.env["CXX"]])
    if "CC" not in config.env.arguments: config.env["CC"] = "mpicc.mpich"
    if "LINK" not in config.env.arguments: config.env["LINK"] = "mpicxx.mpich"
    if "CXX" not in config.env.arguments: config.env["CXX"] = "mpicxx.mpich"
    Helper.printInfo("Using mpich.")

  versionString = getOutput([config.env["CXX"], "-dumpversion"])
  if "." not in versionString:
    versionString = getOutput([config.env["CXX"], "-dumpfullversion"])
  version = config.env._get_major_minor_revision(versionString)
  Helper.printInfo("Using {} {}".format(config.env["CXX"], versionString))

  if not config.CheckExec(config.env["CXX"]) or not config.CheckExec(config.env["CC"]) or \
      not config.CheckExec(config.env["LINK"]) :
    Helper.printErrorAndExit("Compiler executable not found!")

  if not config.CheckCompiler():
    Helper.printErrorAndExit("Compiler found, but it is not working! (Hint: check flags)")

  allWarnings = \
      "-Wall -Wextra \
      -Wcast-qual -Wconversion -Wformat=2 \
      -Wformat-nonliteral -Wformat-security -Winit-self  \
      -Wmissing-format-attribute \
      -Wmissing-include-dirs -Wpacked \
      -Wunreachable-code -Wunused \
      -Wno-unused-parameter".split(" ")

  if not config.env['USE_HPX']:
    allWarnings.append(['-Wswitch-enum', '-Wredundant-decls', '-pedantic'])
  else:
    allWarnings.append(['-Wno-conversion', '-Wno-format-nonliteral'])


  # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or
  #     http://www.swig.org/Release/CHANGES, 03/02/2006
  #     "If you are going to use optimizations turned on with gcc > 4.0 (for example -O2),
  #     ensure you also compile with -fno-strict-aliasing"
  config.env.Append(CPPFLAGS=allWarnings + [
      "-fno-strict-aliasing",
      "-funroll-loops", "-mfpmath=sse"])

  # Mitigation for old Ubuntu (should probably be also applied to Debian?):
  # Package 'libomp-dev' installs a symlink 'libgomp.so' to 'libomp.so' in /usr/lib/x86_64-linux.
  # If this path is manually added (-L...), then ld uses this symlink and, thus, links against
  # the wrong OpenMP library.
  # The mitigation is to ask gcc for its LIBRARY_PATH in combination with -fopenmp and manually
  # add this as the very first LIBPATH.
  # Note, that this code also disables OpenMP if the mitigation command did not produce an
  # adequate path.
  ubuntuVersion = None

  try:
    match = re.search(r"Ubuntu ([0-9]{2}\.[0-9]{2})", getOutput(["lsb_release", "-d"]))
    if match is not None: ubuntuVersion = match.group(1)
  except OSError:
    pass

  if ubuntuVersion in ["16.04", "16.10", "17.04", "17.10", "18.04", "18.10"]:
    output = getOutput([config.env["CXX"], "-v", "-fopenmp", "-xc", "/dev/null"])
    match = re.search(r"LIBRARY_PATH=(.*?):", output)
    firstLibPath = (match.group(1) if match is not None else None)
    if (firstLibPath is None) or (not os.path.exists(firstLibPath)):
      Helper.printWarning("Mitigation for old Ubuntu failed. Did not get libpath. "
                          "Continuing WITHOUT OpenMP.")
    else:
      Helper.printInfo("Mitigation for old Ubuntu: Manually adding {} as first libpath.".format(
          firstLibPath))
      config.env.Append(LIBPATH=[firstLibPath])
      # Safety first: Manually specity libgomp.so.1 as additional library before -fopenmp
      config.env.Append(LINKFLAGS=["-l:libgomp.so.1", "-fopenmp"])
      config.env.Append(CPPFLAGS=["-fopenmp"])
  else:
    config.env.Append(CPPFLAGS=["-fopenmp"])
    config.env.Append(LINKFLAGS=["-fopenmp"])


  #   # limit the number of errors display to something reasonable (useful for templated code)
  #   config.env.Append(CPPFLAGS=["-fmax-errors=5"])

  # required for profiling
  config.env.Append(CPPFLAGS=["-fno-omit-frame-pointer"])

  # GCC has support for colored output since 4.9
  if (version >= (4, 9, 0)) and Helper.terminalSupportsColors():
    config.env.Append(CPPFLAGS=["-fdiagnostics-color=always"])

  if config.env["BUILD_STATICLIB"]:
    config.env.Append(CPPFLAGS=["-D_BUILD_STATICLIB"])

  if config.env["ARCH"] == "sse3":
    config.env.AppendUnique(CPPFLAGS=["-msse3"])
  elif config.env["ARCH"] == "sse42":
    config.env.AppendUnique(CPPFLAGS=["-msse4.2"])
  elif config.env["ARCH"] == "avx":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
  elif config.env["ARCH"] == "fma4":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
    config.env.AppendUnique(CPPFLAGS=["-mfma4"])
  elif config.env["ARCH"] == "avx2":
    config.env.AppendUnique(CPPFLAGS=["-mavx2"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  elif config.env["ARCH"] == "avx512":
    config.env.AppendUnique(CPPFLAGS=["-mavx512f"])
    config.env.AppendUnique(CPPFLAGS=["-mavx512cd"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  else:
    Helper.printErrorAndExit("You must specify a valid ARCH value for gnu.",
                             "Available configurations are: sse3, sse42, avx, fma4, avx2, avx512")

  # check if using MinGW (g++ on win32)
  if config.env["PLATFORM"] == "win32":
    # disable warnings which occur when including Boost in the tests
    # note that definition of hypot is necessary for to the current version of
    # mingw (6.3) and the python interface (see http://stackoverflow.com/questions/10660524/error-building-boost-1-49-0-with-gcc-4-7-0)
    # -> could be removed in the future hopefully
    config.env.Append(CPPFLAGS=["-Wno-switch-enum", "-Wno-deprecated-declarations", "-D_hypot=hypot"])
    # also use "lib" prefix on MinGW for consistency with Linux (default is no prefix)
    config.env["SHLIBPREFIX"] = "lib"

def configureClangCompiler(config):
  if "CC" not in config.env.arguments: config.env["CC"] = "clang"
  if "LINK" not in config.env.arguments: config.env["LINK"] = "clang++"
  if "CXX" not in config.env.arguments: config.env["CXX"] = "clang++"

  versionString = getOutput([config.env["CXX"], "--version"])
  try:
    versionString = re.match(r"^.*version ([^ ]+)", versionString).group(1)
  except AttributeError:
    versionString = "(unknown version)"
  Helper.printInfo("Using {} {}".format(config.env["CXX"], versionString))

  if not config.CheckExec(config.env["CXX"]) or not config.CheckExec(config.env["CC"]) or \
      not config.CheckExec(config.env["LINK"]) :
    Helper.printErrorAndExit("Compiler not found!")

  if not config.CheckCompiler():
    Helper.printErrorAndExit("Compiler found, but it is not working! (Hint: check flags)")

  allWarnings = [
    "-Wall", "-Wextra", "-Weverything", "-Wno-unknown-warning-option",
    "-Wno-c++98-compat-local-type-template-args", "-Wno-c++98-compat-pedantic", "-Wno-deprecated",
    "-Wno-disabled-macro-expansion", "-Wno-documentation", "-Wno-documentation-deprecated-sync",
    "-Wno-documentation-unknown-command", "-Wno-exit-time-destructors", "-Wno-float-equal",
    "-Wno-global-constructors", "-Wno-missing-noreturn", "-Wno-missing-prototypes", "-Wno-padded",
    "-Wno-shadow", "-Wno-shadow-field", "-Wno-sign-conversion", "-Wno-undef",
    "-Wno-unused-parameter", "-Wno-weak-vtables", "-Wno-used-but-marked-unused",
  ]
  config.env.Append(CPPFLAGS=allWarnings)

  config.env.Append(CPPFLAGS=["-fopenmp"])
  config.env.Append(LINKFLAGS=["-fopenmp"])

  if config.env["BUILD_STATICLIB"]:
    config.env.Append(CPPFLAGS=["-D_BUILD_STATICLIB"])

  if config.env["ARCH"] == "sse3":
    config.env.AppendUnique(CPPFLAGS=["-msse3"])
  elif config.env["ARCH"] == "sse42":
    config.env.AppendUnique(CPPFLAGS=["-msse4.2"])
  elif config.env["ARCH"] == "avx":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
  elif config.env["ARCH"] == "fma4":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
    config.env.AppendUnique(CPPFLAGS=["-mfma4"])
  elif config.env["ARCH"] == "avx2":
    config.env.AppendUnique(CPPFLAGS=["-mavx2"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  elif config.env["ARCH"] == "avx512":
    config.env.AppendUnique(CPPFLAGS=["-mavx512f"])
    config.env.AppendUnique(CPPFLAGS=["-mavx512cd"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  else:
    Helper.printErrorAndExit("You must specify a valid ARCH value for clang.",
                             "Available configurations are: sse3, sse4.2, avx, fma4, avx2, avx512")

def configureIntelCompiler(config):
  config.env.AppendUnique(CPPFLAGS=["-Wall", "-ansi", "-Wno-deprecated", "-wd1125",
                                    "-fno-strict-aliasing",
                                    "-ip", "-ipo", "-funroll-loops",
                                    "-ansi-alias", "-fp-speculation=safe",
                                    "-qno-offload"])
  if config.env["COMPILER"] == "intel.mpi":
    if "CC" not in config.env.arguments: config.env["CC"] = "mpiicc"
    if "LINK" not in config.env.arguments: config.env["LINK"] = "mpiicpc"
    if "CXX" not in config.env.arguments: config.env["CXX"] = "mpiicpc"
    config.env["CPPDEFINES"]["USE_MPI"] = "1"
  else:
    if "CC" not in config.env.arguments: config.env["CC"] = "icc"
    if "LINK" not in config.env.arguments: config.env["LINK"] = "icpc"
    if "CXX" not in config.env.arguments: config.env["CXX"] = "icpc"

  versionString = getOutput([config.env["CXX"], "-dumpversion"])
  Helper.printInfo("Using {} {}".format(config.env["CXX"], versionString))

  if not config.CheckExec(config.env["CXX"]) or not config.CheckExec(config.env["CC"]) or \
      not config.CheckExec(config.env["LINK"]) :
    Helper.printErrorAndExit("Compiler not found!")

  if not config.CheckCompiler():
    Helper.printErrorAndExit("Compiler found, but it is not working! (Hint: check flags)")

#   if not config.env["USE_HPX"]:
  config.env.AppendUnique(CPPFLAGS=["-qopenmp"])
  config.env.AppendUnique(LINKFLAGS=["-qopenmp"])

  if config.env["BUILD_STATICLIB"]:
    config.env.AppendUnique(CPPFLAGS=["-D_BUILD_STATICLIB"])

  if config.env["ARCH"] == "sse3":
    config.env.AppendUnique(CPPFLAGS=["-msse3"])
  elif config.env["ARCH"] == "sse42":
    config.env.AppendUnique(CPPFLAGS=["-msse4.2"])
  elif config.env["ARCH"] == "avx":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
  elif config.env["ARCH"] == "avx2":
    config.env.AppendUnique(CPPFLAGS=["-xCORE-AVX2"])
    config.env.AppendUnique(CPPFLAGS=["-fma"])
  elif config.env["ARCH"] == "avx512":
    config.env.AppendUnique(CPPFLAGS=["-xCOMMON-AVX512"])
    config.env.AppendUnique(CPPFLAGS=["-fma"])
  elif config.env["ARCH"] == "mic":
    config.env.AppendUnique(CPPFLAGS=["-mmic"])
    config.env.AppendUnique(LINKFLAGS=["-mmic"])
    config.env["CPPDEFINES"]["USEMIC"] = "1"
  else:
    Helper.printErrorAndExit("You must specify a valid ARCH value for intel.",
                             "Available configurations are: sse3, sse4.2, avx, avx2, avx512, mic")

def detectGSL(config):
  if "GSL_INCLUDE_PATH" in config.env:
    config.env.AppendUnique(CPPPATH=[config.env["GSL_INCLUDE_PATH"]])
  if "GSL_LIBRARY_PATH" in config.env:
    config.env.AppendUnique(LIBPATH=[config.env["GSL_LIBRARY_PATH"]])
  if config.CheckCXXHeader("gsl/gsl_version.h") and config.CheckLib(["gsl", "gslcblas"], language="c++", autoadd=0):
    Helper.printInfo("GSL is installed, enabling GSL support.")
    config.env["USE_GSL"] = True
    config.env["CPPDEFINES"]["USE_GSL"] = "1"
  elif config.env["USE_GSL"]:
    Helper.printErrorAndExit("gsl/gsl_version.h or libsgl/libgslcblas were not found, but required for GSL")
  else:
    Helper.printInfo("GSL support could not be enabled.")

def detectZlib(config):
  if config.CheckLib("z", language="c++", autoadd=0) and config.CheckCXXHeader("zlib.h"):
    Helper.printInfo("zlib is installed, enabling ZLIB support.")
    config.env["USE_ZLIB"] = True
    config.env["CPPDEFINES"]["ZLIB"] = "1"
  elif config.env["USE_ZLIB"]:
    if config.env["PLATFORM"] == "win32":
      Helper.printWarning("zlib is currently not supported on Windows. Continuing withouth zlib.")
    else:
      Helper.printErrorAndExit("USE_ZLIB is set but either libz or zlib.h is missing!")
  else:
    Helper.printInfo("ZLIB support could not be enabled.")

def detectScaLAPACK(config):
  if "SCALAPACK_LIBRARY_PATH" in config.env:
    config.env.AppendUnique(LIBPATH=[config.env["SCALAPACK_LIBRARY_PATH"]])

  # check if USE_SCALAPACK was given as parameter and is disabled
  if "USE_SCALAPACK" in config.env and not config.env["USE_SCALAPACK"]:
    Helper.printInfo("ScaLAPACK disabled.")
    return

  if config.env["COMPILER"] not in ("openmpi", "mpich", "intel.mpi"):
    if "USE_SCALAPACK" in config.env and config.env["USE_SCALAPACK"]:
      Helper.printErrorAndExit("USE_SCALAPACK was set, but no mpi compiler was used (openmpi, mpich or intel.mpi)")
    config.env["USE_SCALAPACK"] = False
    Helper.printInfo("ScaLAPACK was disabled as no mpi compiler was used (openmpi, mpich or intel.mpi)")
    return


  # check if there is a ScaLAPACK version installed
  if "SCALAPACK_LIBRARY_NAME" in config.env and config.CheckLib(config.env["SCALAPACK_LIBRARY_NAME"], language="c++", autoadd=0):
    config.env["SCALAPACK_VERSION"] = "custom"
    config.env["USE_SCALAPACK"] = True
    config.env["CPPDEFINES"]["USE_SCALAPACK"] = "1"
    Helper.printInfo("Using scalapack version from SCALAPACK_LIBRARY_NAME: " + str(config.env["SCALAPACK_LIBRARY_NAME"]))
  elif config.CheckMKL():
    config.env["SCALAPACK_VERSION"] = "mkl"
    config.env["USE_SCALAPACK"] = True
    config.env["CPPDEFINES"]["USE_SCALAPACK"] = "1"
    Helper.printInfo("Using mkl ScaLAPACK")
  elif config.CheckLib("scalapack", language="c++", autoadd=0):
    config.env["SCALAPACK_VERSION"] = "netlib"
    config.env["USE_SCALAPACK"] = True
    config.env["CPPDEFINES"]["USE_SCALAPACK"] = "1"
    Helper.printInfo("Using netlib ScaLAPACK")
  elif config.env["COMPILER"] == "openmpi" and config.CheckLib("scalapack-openmpi", language="c++", autoadd=0):
    config.env["SCALAPACK_VERSION"] = "openmpi"
    config.env["USE_SCALAPACK"] = True
    config.env["CPPDEFINES"]["USE_SCALAPACK"] = "1"
    Helper.printInfo("Using openmpi ScaLAPACK")
  elif config.env["COMPILER"] == "mpich" and config.CheckLib("scalapack-mpich", language="c++", autoadd=0):
    config.env["SCALAPACK_VERSION"] = "mpich"
    config.env["USE_SCALAPACK"] = True
    config.env["CPPDEFINES"]["USE_SCALAPACK"] = "1"
    Helper.printInfo("Using mpich ScaLAPACK")
  elif "USE_SCALAPACK" in config.env and config.env["USE_SCALAPACK"]:
    Helper.printErrorAndExit("No supported version of ScaLAPACK was found.")
  else:
    config.env["USE_SCALAPACK"] = False
    Helper.printInfo("No ScaLAPACK version found, ScaLAPACK support disabled.")
