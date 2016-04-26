# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from __future__ import print_function

import io
import glob
import os
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
    if module in ["jsgpp", "pysgpp"]:
      languageSupport.append(module)
      continue
    modules.append(module)

  modules.sort()
  languageSupport.sort()
  return modules, languageSupport

# Definition of flags / command line parameters for SCons
#########################################################################

def multiParamConverter(s):
  return s.split(",")

# this class is used with "sys.stdout = Logger(sys.stdout)" at the very beginning
# => print lines to the terminal and to the log file simultaneously
class Logger(object):
  def __init__(self, terminal):
    self.terminal = terminal

    # clear file
    with open("build.log", "a") as logFile:
      logFile.seek(0)
      logFile.truncate()

  def write(self, message):
    self.terminal.write(message)
    # use io.open instead of open because Python replaces all "\n" with os.linesep by default,
    # i.e., on Windows, if we call write with some "\r\n" in it, they get replaced by "\r\r\n";
    # newline="" prevents this conversion
    with io.open("build.log", "a", newline="") as logFile:
      logFile.write(unicode(message))

  def flush(self):
    self.terminal.flush()

class FinalMessagePrinter(object):
  def __init__(self):
    self.enabled = True
    self.env = None
    self.sgppBuildPath = None
    self.pysgppPackagePath = None
    self.exitCode = None
    self.exitFunction = sys.exit
    sys.exit = self.exit

  def enable(self):
    self.enabled = True

  def disable(self):
    self.enabled = False

  def exit(self, exitCode=0):
    self.exitCode = exitCode
    self.exitFunction(exitCode)

  def printMessage(self):
    if not self.enabled:
      return

    failures = SCons.Script.GetBuildFailures()
    build_was_successful = (len(failures) == 0) and (SCons.Script.Main.exit_status == 0)
    build_interrupted = build_was_successful and (self.exitCode != 0)

    if build_interrupted:
      return

    if build_was_successful:
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
      print("""# ------------------------------------------------------------------------
# An error occurred while compiling SG++.
# If you believe this is a bug in SG++, please attach the build.log and
# the config.log file when contacting the SG++ developers.
# ------------------------------------------------------------------------""")

# detour compiler output to print dots if VERBOSE=0
def printCommand(s, target, src, env):
  if env["VERBOSE"]:
    sys.stdout.write(u"%s\n" % s)
  else:
    sys.stdout.write(u".")
    sys.stdout.flush()

#creates a Doxyfile containing proper module paths based on Doxyfile_template
def prepareDoxyfile(modules):
  '''Create Doxyfile(s) and overview-pages
  @param modules list of modules'''

  # create Doxyfile
  with open("Doxyfile_template", "r") as doxyFileTemplate:
    with open("Doxyfile", "w") as doxyFile:
      inputPaths = "INPUT ="
      excludePaths = "EXCLUDE ="
      examplePaths = "EXAMPLE_PATH ="
      imagePaths = "IMAGE_PATH ="

      for moduleName in modules:
        inputPath = moduleName + "/"
        examplePath = moduleName + "/examples"
        testPath = moduleName + "/tests"
        imagePath = moduleName + "/doc/doxygen/images"

        #print(os.path.join(os.getcwd(),inputPath))
        if os.path.exists(os.path.join(os.getcwd(), inputPath)):
          inputPaths += " " + inputPath
        if os.path.exists(os.path.join(os.getcwd(), examplePath)):
          examplePaths += " " + examplePath
          excludePaths += " " + examplePath
        if os.path.exists(os.path.join(os.getcwd(), testPath)):
          excludePaths += " " + testPath
        if os.path.exists(os.path.join(os.getcwd(), imagePath)):
          imagePaths += " " + imagePath

      for line in doxyFileTemplate.readlines():
        if re.match(r"INPUT  .*", line):
          doxyFile.write(inputPaths + "\n")
        elif re.match(r"EXCLUDE  .*", line):
          doxyFile.write(excludePaths + "\n")
        elif re.match(r"EXAMPLE_PATH  .*", line):
          doxyFile.write(examplePaths + "\n")
        elif re.match(r"IMAGE_PATH  .*", line):
          doxyFile.write(imagePaths + "\n")
        else:
          doxyFile.write(line)

  # create example menu page
  with open("base/doc/doxygen/examples.doxy", "w") as examplesFile:
    examplesFile.write('''/**
@page examples Examples

This is a collection of examples from all modules.

If you're new to SG++ or want to try out quickly,
read the @ref code_examples_tutorial first.

To add new examples to the documentation,
go to the respective folder MODULE_NAME/doc/doxygen/ and
add a new example file code_examples_NAME.doxy with doxygen-internal
name code_examples_NAME.

Note that SCons automatically compiles (but not runs)
all C++ examples on each run.
For this to work, the examples must lie in the directories of the form
\c /path/to/SGpp/trunk/MODULE_NAME/examples.

''')

    modules.sort()
    tutorial = "code_examples_tutorial"

    for moduleName in modules:
      examplesFile.write("<h2>Module {}</h2>\n".format(moduleName))
      subpages = glob.glob(os.path.join(
        moduleName, "doc", "doxygen", "code_examples_*.doxy"))
      subpages = [os.path.split(path)[-1][:-5]
            for path in glob.glob(os.path.join(
              moduleName, "doc", "doxygen",
              "code_examples_*.doxy"))]
      subpages.sort()
      if tutorial in subpages:
        del subpages[subpages.index(tutorial)]
        subpages = [tutorial] + subpages

      for subpage in subpages:
        examplesFile.write("- @subpage {}\n".format(subpage))

    examplesFile.write("**/\n")

  # create module page
  with open("base/doc/doxygen/modules.doxy", "w") as modulesFile:
    with open("base/doc/doxygen/modules.stub0", "r") as stubFile:
      modulesFile.write(stubFile.read())

    for moduleName in modules:
      for subpage in glob.glob(os.path.join(moduleName, "doc", "doxygen", "module_*.doxy")):
        modulesFile.write("- @subpage " + os.path.splitext(os.path.split(subpage)[-1])[0] + "\n")

    with open("base/doc/doxygen/modules.stub1", "r") as stubFile:
      modulesFile.write(stubFile.read())


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
      for var in spawnEnv:
        spawnEnv[var] = spawnEnv[var].encode("ascii", "replace")

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
  import msvcrt
  import win32api
  import win32con
  import win32event
  import win32file
  import win32pipe
  import win32process
  import win32security

  def win32Spawn(sh, escape, cmd, args, spawnEnv):
    for var in spawnEnv:
      spawnEnv[var] = spawnEnv[var].encode("ascii", "replace")

    sAttrs = win32security.SECURITY_ATTRIBUTES()
    sAttrs.bInheritHandle = 1
    hStdinR, hStdinW = win32pipe.CreatePipe(sAttrs, 0)
    hStdoutR, hStdoutW = win32pipe.CreatePipe(sAttrs, 0)
    hStderrR, hStderrW = win32pipe.CreatePipe(sAttrs, 0)

    pid = win32api.GetCurrentProcess()

    def replaceHandle(handle) :
      tmp = win32api.DuplicateHandle(pid, handle, pid, 0, 0, win32con.DUPLICATE_SAME_ACCESS)
      win32file.CloseHandle(handle)
      return tmp

    hStdinW = replaceHandle(hStdinW)
    hStdoutR = replaceHandle(hStdoutR)
    hStderrR = replaceHandle(hStderrR)

    sAttrs = win32security.SECURITY_ATTRIBUTES()
    startupInfo = win32process.STARTUPINFO()
    startupInfo.hStdInput = hStdinR
    startupInfo.hStdOutput = hStdoutW
    startupInfo.hStdError = hStderrW
    startupInfo.dwFlags = win32process.STARTF_USESTDHANDLES
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

      win32file.CloseHandle(hStdinR)
      win32file.CloseHandle(hStdoutW)
      win32file.CloseHandle(hStderrW)
      with os.fdopen(msvcrt.open_osfhandle(hStdoutR, 0), "rb") as f: stdout = f.read()
      with os.fdopen(msvcrt.open_osfhandle(hStderrR, 0), "rb") as f: stderr = f.read()
      sys.stdout.write(stdout)
      sys.stderr.write(stderr)

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

