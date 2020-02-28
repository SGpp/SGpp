# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from __future__ import print_function
import glob
import os
import platform
import re
import string
import subprocess
import sys

import SCons.Script
import SCons.Script.Main

def printInfo(*s):
  print("Info:", "\n      ".join([str(x) for x in s]))

def printWarning(*s):
  print("Warning:", "\n         ".join([str(x) for x in s]), file=sys.stderr)

def printErrorAndExit(*s):
  print("\nError:", "\n       ".join([str(x) for x in s]), file=sys.stderr)
  sys.exit(1)

# get all folders containing an "SConscript*" file
# path has to end with "/"
def getModules(ignoreFolders):

#     if path[-1] != "/":
#         path += "/"
  #path = "#/"
  path = ""
  suffix = "/SConscript"
  searchString = path + "*" + suffix
  modulePaths = glob.glob(searchString)
  modules = []
  languageSupport = []

  for modulePath in modulePaths:
    module = modulePath[:-len(suffix)]
    module = module[len(path):]
    if module in ignoreFolders:
      continue
    if module in ["jsgpp", "pysgpp", "matsgpp"]:
      languageSupport.append(module)
      continue
    modules.append(module)

  modules.sort()
  languageSupport.sort()
  return modules, languageSupport

# Definition of flags / command line parameters for SCons
#########################################################################

def multiParamConverter(s):
  return s.split(" ")

def multiParamPathConverter(s):
  return s.split(os.pathsep)

def multiParamDefineConverter(s):
  return (dict([x.split("=") for x in s.split(" ")]) if s != "" else {})

# this class is used with "sys.stdout = Logger(sys.stdout)" at the very beginning
# => print lines to the terminal and to the log file simultaneously
class Logger(object):
  def __init__(self, terminal):
    # save original terminal
    self.terminal = terminal
    # regex for ANSI color codes
    self.ansiCode = re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]")

    # clear file
    with open("build.log", "a") as logFile:
      logFile.seek(0)
      logFile.truncate()

  def write(self, message):
    # write on terminal
    self.terminal.write(message)

    # remove ANSI color codes (e.g., for GCC >= 4.9)
    # (from https://stackoverflow.com/a/33925425)
    message = self.ansiCode.sub("", message)

    # Python replaces all "\n" with os.linesep by default,
    # i.e., on Windows, if we call write with some "\r\n" in it, they get replaced by "\r\r\n"
    message = message.replace(os.linesep, "\n")

    with open("build.log", "a") as logFile:
      logFile.write(message)

  def flush(self):
    self.terminal.flush()

  def isatty(self):
    return (self.terminal.isatty() if hasattr(self.terminal, "isatty") else None)

# class for printing a final success/error message,
# depending on whether there were compilation errors or not
class FinalMessagePrinter(object):
  def __init__(self):
    # whether the final message will be printed or not
    self.enabled = True
    # global SCons environment (has to be set from outside)
    self.env = None
    # path of the SG++ libraries (has to be set from outside)
    self.sgppBuildPath = None
    # path of the pysgpp package (has to be set from outside)
    self.pysgppPackagePath = None
    # exit code if exited with sys.exit
    self.exitCode = None
    # override sys.exit (called by SCons) to save the exit code,
    # since atexit.register doesn't give access to the exit code
    self.exitFunction = sys.exit
    sys.exit = self.exit

  def enable(self):
    # enable printing of final message
    self.enabled = True

  def disable(self):
    # disable printing of final message
    self.enabled = False

  def exit(self, exitCode=0):
    # save exit code, call original sys.exit
    self.exitCode = exitCode
    self.exitFunction(exitCode)

  def printMessage(self):
    # check if we should print
    if not self.enabled:
      return

    # list of SCons build failures
    failures = SCons.Script.GetBuildFailures()
    # the build was successful, if there were no failures and the exit_status is zero
    build_was_successful = (len(failures) == 0) and (SCons.Script.Main.exit_status == 0)
    # the build was interrupted, if there were not errors, but the exit code is non-zero
    build_interrupted = build_was_successful and (self.exitCode != 0)

    # don't print if Ctrl+C was hit
    if build_interrupted:
      return

    if build_was_successful:
      # print success message (instructions)
      if self.env["PLATFORM"] in ["cygwin", "win32"]:
        filename = "INSTRUCTIONS_WINDOWS"
      elif self.env["PLATFORM"] == "darwin":
        filename = "INSTRUCTIONS_MAC"
      else:
        filename = "INSTRUCTIONS"

      with open(filename) as f:
        instructionsTemplate = string.Template(f.read())
        print()
        print(instructionsTemplate.safe_substitute(SGPP_BUILD_PATH=self.sgppBuildPath,
                                                   PYSGPP_PACKAGE_PATH=self.pysgppPackagePath))
    else:
      # print error message
      print("""# ------------------------------------------------------------------------
# An error occurred while compiling SG++.
# If you believe this is a bug in SG++, please attach the build.log and
# the config.log file when contacting the SG++ developers.
# ------------------------------------------------------------------------""")

# returns whether the terminal supports ANSI color codes
# (actually, we only check if it's a tty and we're not on Windows,
# since there doesn't seem to be a portable way to check that)
def terminalSupportsColors():
  return hasattr(sys.stderr, "isatty") and sys.stderr.isatty() and (platform.system() != "Windows")

# detour compiler output to print dots if VERBOSE=0
def printCommand(s, target, src, env):
  if env["VERBOSE"]:
    sys.stdout.write(u"%s\n" % s)
  else:
    sys.stdout.write(u".")
    sys.stdout.flush()

def flatDependencyGraph(dependencies, acc):
  for dependency in dependencies[::-1]:
    if dependency not in acc:
      acc = [dependency] + acc
  return acc

# Override env["SPAWN"] to log output (stdout/stderr) of child processes.
# (see https://bitbucket.org/scons/scons/wiki/BuildLog)
def setSpawn(env):
  def echoSpawn(sh, escape, cmd, args, spawnEnv):
      """Spawn which echos stdout/stderr from the child."""
      # convert spawnEnv from unicode strings
      # for var in spawnEnv:
      #   spawnEnv[var] = spawnEnv[var].encode("ascii", "replace")

      newArgs = " ".join(args[1:])
      cmdLine = cmd + " " + newArgs

      p = subprocess.Popen(
          cmdLine,
          env=spawnEnv,
          stderr=subprocess.PIPE,
          stdout=subprocess.PIPE,
          shell=True,
          universal_newlines=True)
      (stdout, stderr) = p.communicate()

      # Does this screw up the relative order of the two?
      sys.stdout.write(stdout)
      sys.stderr.write(stderr)
      return p.returncode

  env["SPAWN"] = echoSpawn

# On win32, the command lines are limited to a ridiculously short length
# (1000 chars). However, compiler/linker command lines easily exceed that
# length. The following is a fix for that.
# It has to be enabled with "env["SPAWN"] = win32Spawn".
# (see https://bitbucket.org/scons/scons/wiki/LongCmdLinesOnWin32)
def setWin32Spawn(env):
  # the lines commented out were to write the results of the child processes to build.log
  # (not only the commands which were executed)
  # this seems, however, to cause sometimes deadlocks when compiling with -j
  #import msvcrt
  import win32api
  #import win32con
  import win32event
  import win32file
  #import win32pipe
  import win32process
  import win32security

  def win32Spawn(sh, escape, cmd, args, spawnEnv):
    #sAttrs = win32security.SECURITY_ATTRIBUTES()
    #sAttrs.bInheritHandle = 1
    #hStdinR, hStdinW = win32pipe.CreatePipe(sAttrs, 0)
    #hStdoutR, hStdoutW = win32pipe.CreatePipe(sAttrs, 0)
    #hStderrR, hStderrW = win32pipe.CreatePipe(sAttrs, 0)
    #
    #pid = win32api.GetCurrentProcess()
    #
    #def replaceHandle(handle) :
    #  tmp = win32api.DuplicateHandle(pid, handle, pid, 0, 0, win32con.DUPLICATE_SAME_ACCESS)
    #  win32file.CloseHandle(handle)
    #  return tmp
    #
    #hStdinW = replaceHandle(hStdinW)
    #hStdoutR = replaceHandle(hStdoutR)
    #hStderrR = replaceHandle(hStderrR)

    sAttrs = win32security.SECURITY_ATTRIBUTES()
    startupInfo = win32process.STARTUPINFO()
    #startupInfo.hStdInput = hStdinR
    #startupInfo.hStdOutput = hStdoutW
    #startupInfo.hStdError = hStderrW
    #startupInfo.dwFlags = win32process.STARTF_USESTDHANDLES
    newArgs = " ".join(map(escape, args[1:]))
    cmdLine = cmd + " " + newArgs

    # check for any special operating system commands
    if cmd == "del":
      for arg in args[1:]:
        win32file.DeleteFile(arg)
      exitCode = 0
    else:
      # otherwise execute the command.
      try:
        hProcess, hThread, dwPid, dwTid = win32process.CreateProcess(
            None, cmdLine, None, None, 1, 0, spawnEnv, None, startupInfo)
      except:
        errorCode = win32api.GetLastError()
        raise RuntimeError("Could not execute the following " +
          "command line (error code {}): {}".format(errorCode, cmdLine))
      win32event.WaitForSingleObject(hProcess, win32event.INFINITE)
      exitCode = win32process.GetExitCodeProcess(hProcess)

      #win32file.CloseHandle(hStdinR)
      #win32file.CloseHandle(hStdoutW)
      #win32file.CloseHandle(hStderrW)
      #with os.fdopen(msvcrt.open_osfhandle(hStdoutR, 0), "rb") as f: stdout = f.read()
      #with os.fdopen(msvcrt.open_osfhandle(hStderrR, 0), "rb") as f: stderr = f.read()
      #sys.stdout.write(stdout)
      #sys.stderr.write(stderr)

      win32file.CloseHandle(hProcess)
      win32file.CloseHandle(hThread)
    return exitCode

  env["SPAWN"] = win32Spawn

# get all subdirs of path, required by CheckJNI
def getSubdirs(path):
  pathList = []
  for f in os.listdir(path):
    if os.path.isdir(os.path.join(path, f)):
      pathList.append(os.path.join(path, f))
  return pathList

# Custom test for executables used during configuration
def CheckExec(context, cmd):
  context.Message("Checking for %s..." % (cmd))
  ret = context.env.WhereIs(cmd)
  if ret == None:
    ret = ""
  context.Result(ret)
  return ret

def CheckCompiler(context):

  checkSrcFile = """
int main(int argc, char **argv) {
  /**
   * Does nothing, just a test whether the compile works.
   */
}
  """
  context.Message("Testing whether compiler works... ")
  result = context.TryCompile(checkSrcFile, ".cpp")
  context.Result(result)
  return result

def CheckFlag(context, flagString):

  checkSrcFile = """
int main(int argc, char **argv) {
  /**
   * Does nothing, just a test whether a compilation with a specific flag is possible.
   */
}
  """
  context.Message("Checking for flag \"" + flagString + "\"... ")
  oldFlags = context.env["CPPFLAGS"]
  context.env.AppendUnique(CPPFLAGS = flagString)
  result = context.TryCompile(checkSrcFile, ".cpp")
  context.Result(result)
  context.env.Replace(CPPFLAGS = oldFlags)
  return result

# Check for jni header file
# if found, additionally add all subdirs to CPPPATH (platform dependent files)
def CheckJNI(context):
  print("Trying to locate jni.h...")
  # message if JNI_CPPINCLUDE not set
  if not os.environ.get("JNI_CPPINCLUDE"):
    print("... JNI_CPPINCLUDE not set")
  # check for JAVA_HOME first
  if os.environ.get("JAVA_HOME"):
    pname = os.path.join(os.environ.get("JAVA_HOME"), "include")
    if os.path.exists(os.path.join(pname, "jni.h")):
      context.env.Append(CPPPATH=[pname] + getSubdirs(pname))
      res = "... found in " + pname
      context.Result(res)
      return res
    else:
      print("... not found in $JAVA_HOME/include")
  else:
    print("... JAVA_HOME not set")

  # not found, try guessing:
  # look, where java and javac are located:
  # include/ directory might be 1 or 2 dirs below
  print("... trying to guess")
  for f in ["java", "javacc"]:
    fdir = context.env.WhereIs(f)
    if not fdir:
      continue
    # os.path.realpath to resolve links
    baseDir = os.path.dirname(os.path.realpath(fdir))
    for subDir in ["..", os.path.join("..", "..")]:
      pname = os.path.join(baseDir, subDir, "include")
      if os.path.exists(os.path.join(pname, "jni.h")):
        context.env.Append(CPPPATH=[pname] + getSubdirs(pname))
        res = "... found in " + pname
        context.Result(res)
        return res

  context.Result("... nothing found!")
  return 0

# check for the mkl scalapack version
# (Conftest CheckLib version with extra_libs argument is needed)
def CheckMklScalapack(context):
  mkl_libs = ["mkl_blacs_intelmpi_lp64",
                "mkl_gf_lp64",
                "mkl_gnu_thread",
                "mkl_core",
                "gomp",
                "pthread",
                "dl",
                "m"]
  res = not SCons.Conftest.CheckLib(
    context, ["mkl_scalapack_lp64"], extra_libs=mkl_libs, language="c++", autoadd=0)
  context.did_show_result = 1
  context.Result(res)
  return res

# check if library is installed and set define flag,
# use None for "libraries" if it's a header-only library,
# "additionalDependencies" is a list which will be extended by the libraries
def checkForLibrary(config, env, additionalDependencies,
                    name, flag, headers, libraries):
  if type(headers) is str: headers = [headers]
  if type(libraries) is str: libraries = [libraries]
  if libraries is None: libraries = []

  if (flag not in config.env) or (not config.env[flag]):
    printInfo("SG++ will be compiled without {} (flag not set).".format(name))
  elif not config.CheckHeader(headers, language="C++"):
    printErrorAndExit("The flag {} was given, but the".format(flag),
                      "necessary headers {} were not found.".format(headers))
  elif (len(libraries) > 0) and (not config.CheckLib(libraries, language="C++")):
    printErrorAndExit("The flag {} was given, but the".format(flag),
                      "necessary libraries {} were not found.".format(libraries))
  else:
    printInfo("SG++ will be compiled with {}.".format(name))
    additionalDependencies.extend(libraries)
    env["CPPDEFINES"][flag] = "1"
