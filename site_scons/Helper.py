# Copyright (C) 2008-today The SG++ Project
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

def convertExampleSourceToDoxy(sourcePath):
  with open(sourcePath, "r") as f: source = f.read()

  sourcePath = os.path.normpath(sourcePath)
  sourcePathComponents = sourcePath.split(os.sep)
  sourceFileName = sourcePathComponents[-1]
  sourceFileType = os.path.splitext(sourceFileName)[1][1:]
  moduleName = sourcePathComponents[sourcePathComponents.index("examples") - 1]

  if sourceFileType == "cpp":
    doxygenBlockCommentBegin = "/**"
    doxygenBlockCommentEnd = "*/"
    doxygenLineCommentBegin = "///"
  elif sourceFileType == "py":
    doxygenBlockCommentBegin = None
    doxygenBlockCommentEnd = None
    doxygenLineCommentBegin = "##"
  elif sourceFileType == "java":
    if " main(" in source:
      doxygenBlockCommentBegin = "/**"
      doxygenBlockCommentEnd = "*/"
      doxygenLineCommentBegin = "///"
    else:
      return None
  elif sourceFileType == "m":
    doxygenBlockCommentBegin = None
    doxygenBlockCommentEnd = None
    doxygenLineCommentBegin = "%%"

  pageMatch = re.search(r"[@\\]page +([^ ]+) +(.*)$", source, flags=re.MULTILINE)
  doxy = "/**\n"

  if pageMatch is not None:
    pageName = pageMatch.group(1)
    doxy += "\\dontinclude {}\n".format(sourceFileName)
    inDoxygenComment = False
    previousNonBlankLines = []
    skipUntilNextNonBlankLine = False
    sawFirstDoxygenComment = False

    for sourceLine in source.splitlines():
      if sourceLine.strip() == "":
        continue

      if sourceLine.strip() == doxygenBlockCommentBegin:
        inDoxygenComment = True
        if sawFirstDoxygenComment:
          if len(previousNonBlankLines) > 0:
            doxy += "\\until {}\n\n".format(previousNonBlankLines[-1])
          else:
            doxy += "\n"
        else:
          # assume up to now, there was only the copyright notice and no code
          # ==> don't include now
          sawFirstDoxygenComment = True
        previousNonBlankLines = []
        skipUntilNextNonBlankLine = False
      elif (sourceLine.strip() == doxygenBlockCommentEnd) and inDoxygenComment:
        inDoxygenComment = False
        skipUntilNextNonBlankLine = True
      elif sourceLine.strip().startswith(doxygenLineCommentBegin):
        if sawFirstDoxygenComment:
          if len(previousNonBlankLines) > 0:
            doxy += "\\until {}\n\n".format(previousNonBlankLines[-1])
        else:
          sawFirstDoxygenComment = True
        i = sourceLine.index(doxygenLineCommentBegin) + len(doxygenLineCommentBegin)
        doxy += sourceLine[i:] + "\n"
        previousNonBlankLines = []
        skipUntilNextNonBlankLine = True
      else:
        if inDoxygenComment:
          doxy += "{}\n".format(sourceLine.strip(" *"))
        else:
          if skipUntilNextNonBlankLine:
            doxy += "\\skip {}\n".format(sourceLine)
            skipUntilNextNonBlankLine = False
          previousNonBlankLines.append(sourceLine)

    if len(previousNonBlankLines) > 0:
      doxy += "\\until {}\n\n".format(previousNonBlankLines[-1])
  else:
    pageName = "example_{}".format(sourceFileName.replace(".", "_"))
    doxy += "@page {} {}\n".format(pageName, sourceFileName)
    doxy += "This example can be found under <tt>{}</tt>.\n".format(sourcePath)
    doxy += "\\include {}\n".format(sourceFileName)

  doxy += "*/\n"
  doxyPath = "{}/doc/doxygen/{}.doxy".format(moduleName, pageName)

  with open(doxyPath, "w") as f: f.write(doxy)
  return {"pageName" : pageName, "language" : sourceFileType, "moduleName" : moduleName}

#creates a Doxyfile containing proper module paths based on Doxyfile_template
def prepareDoxyfile(modules):
  '''Create Doxyfile(s) and overview-pages
  @param modules list of modules'''

  examples = []

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

        for exampleFileName in os.listdir(examplePath):
          if any([exampleFileName.endswith(ext) for ext in [".cpp", ".py", ".java", ".m"]]):
            example = convertExampleSourceToDoxy("{}/{}".format(examplePath, exampleFileName))
            if example is not None:
              examples.append(example)

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

  for language in ["cpp", "py", "java", "m"]:
    examplesInLanguage = [example for example in examples if example["language"] == language]

    # create examples menu page
    with open("base/doc/doxygen/examples_{}.doxy".format(language), "w") as examplesFile:
      examplesFile.write("/**\n")
      languageName = {"cpp" : "C++", "py" : "Python", "java" : "Java", "m" : "MATLAB"}[language]
      examplesFile.write("@page examples_{} {} Examples\n".format(language, languageName))
      examplesFile.write("This is a list of all {} examples.\n".format(languageName))
      examplesFile.write("If you don't know where to start, look at the @ref " +
                         "example_tutorial_{} example first.\n".format(language))
      examplesFile.write("All examples can be found in the <tt>MODULE_NAME/example/</tt> " +
                         "directories.\n")
      if language == "cpp":
        examplesFile.write("Note that SCons automatically compiles (but not runs) " +
                           "all C++ examples on each run. The executables can be found in the " +
                           "same directory in which the examples reside and can be run " +
                           "directly, if <tt>LD_LIBRARY_PATH</tt> (on Linux/Mac) or " +
                           "<tt>PATH</tt> (on Windows) is set correctly.\n")
      examplesFile.write("\nFor more instructions on how to run the examples, " +
                         "please see @ref installation.\n\n")

      moduleNames = sorted(list(set([example["moduleName"] for example in examplesInLanguage])))
      for moduleName in moduleNames:
        examplesInLanguageAndModule = [example for example in examplesInLanguage
                                       if example["moduleName"] == moduleName]
        examplesInLanguageAndModule.sort(key=lambda example: example["pageName"].lower())
        examplesFile.write("\n<b>Module sgpp::{}</b>\n\n".format(moduleName))
        for example in examplesInLanguageAndModule:
          examplesFile.write("- @subpage {}\n".format(example["pageName"]))

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
    for var in spawnEnv:
      spawnEnv[var] = spawnEnv[var].encode("ascii", "replace")

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

