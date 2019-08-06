# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


import glob
import os
import re

import Helper

# convert a example source file to a *.doxy documentation
def convertExampleSourceToDoxy(sourcePath):
  # read source
  with open(sourcePath, "r") as f: source = f.read()

  sourcePath = os.path.normpath(sourcePath)
  sourcePathComponents = sourcePath.split(os.sep)
  sourceFileName = sourcePathComponents[-1]
  sourceFileType = os.path.splitext(sourceFileName)[1][1:]
  moduleName = sourcePathComponents[sourcePathComponents.index("examples") - 1]
  doxyFolder = os.path.join(moduleName, "doc", "doxygen")

  # feasible languages
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
      # skip file if there is no main method
      return None
  elif sourceFileType == "m":
    doxygenBlockCommentBegin = None
    doxygenBlockCommentEnd = None
    doxygenLineCommentBegin = "%%"

  # look if there is a \page command
  pageMatch = re.search(r"[@\\]page +([^ ]+) +(.*)$", source, flags=re.MULTILINE)
  doxy = "/**\n"

  if pageMatch is not None:
    # page command found ==> assume this is a documented example
    pageName = pageMatch.group(1)
    inDoxygenComment = False
    snippet = []
    skipUntilNextNonBlankLine = False
    sawFirstDoxygenComment = False
    snippetNumber = 1

    def saveSnippet(snippetNumber):
      snippetName = "{}_snippet{}.{}".format(
          pageName, snippetNumber, sourceFileType)
      snippetPath = os.path.join(doxyFolder, "snippets", snippetName)
      with open(snippetPath, "w") as f: f.write("\n".join(snippet))
      return snippetName

    for sourceLine in source.splitlines():
      # skip empty lines
      if sourceLine.strip() == "":
        if not skipUntilNextNonBlankLine: snippet.append(sourceLine)
        continue

      if sourceLine.strip() == doxygenBlockCommentBegin:
        # begin of Doxygen block comment
        inDoxygenComment = True
        if sawFirstDoxygenComment:
          if len(snippet) > 0:
            snippetPath = saveSnippet(snippetNumber)
            snippetNumber += 1
            snippet = []
            doxy += "\\include {}\n\n".format(snippetPath)
          else:
            doxy += "\n"
        else:
          # assume up to now, there was only the copyright notice and no code ==> skip
          sawFirstDoxygenComment = True
        snippet = []
        skipUntilNextNonBlankLine = False
      elif (sourceLine.strip() == doxygenBlockCommentEnd) and inDoxygenComment:
        # end of Doxygen block comment
        inDoxygenComment = False
        skipUntilNextNonBlankLine = True
      elif sourceLine.strip().startswith(doxygenLineCommentBegin):
        # Doxygen line comment
        if sawFirstDoxygenComment:
          if len(snippet) > 0:
            snippetPath = saveSnippet(snippetNumber)
            snippetNumber += 1
            snippet = []
            doxy += "\\include {}\n\n".format(snippetPath)
        else:
          sawFirstDoxygenComment = True
        i = sourceLine.index(doxygenLineCommentBegin) + len(doxygenLineCommentBegin)
        doxy += sourceLine[i:] + "\n"
        snippet = []
        skipUntilNextNonBlankLine = True
      elif inDoxygenComment:
        # line of Doxygen block comment
        doxy += "{}\n".format(sourceLine.strip(" *"))
      else:
        # source code line
        skipUntilNextNonBlankLine = False
        snippet.append(sourceLine)

    # write remaining lines if any
    if len(snippet) > 0:
      snippetPath = saveSnippet(snippetNumber)
      snippetNumber += 1
      snippet = []
      doxy += "\\include {}\n\n".format(snippetPath)
  else:
    # no page command ==> include the source verbatim without additional documentation
    pageName = "example_{}".format(sourceFileName.replace(".", "_"))
    doxy += "@page {} {}\n".format(pageName, sourceFileName)
    doxy += "This example can be found under <tt>{}</tt>.\n".format(sourcePath)
    doxy += "\\include {}\n".format(sourceFileName)

  doxy += "*/\n"
  doxyPath = os.path.join(doxyFolder, "{}.doxy".format(pageName))

  # write *.doxy file
  with open(doxyPath, "w") as f: f.write(doxy)
  return {"pageName" : pageName, "language" : sourceFileType, "moduleName" : moduleName}

# search for example source files and convert them to *.doxy files
def convertExampleSourcesToDoxy(modules):
  examples = []

  # for each module
  for moduleName in modules:
    examplePath = os.path.join(moduleName, "examples")
    if not os.path.isdir(examplePath): continue

    snippetsFolder = os.path.join(moduleName, "doc", "doxygen", "snippets")
    if not os.path.isdir(snippetsFolder): os.makedirs(snippetsFolder)

    # search for examples
    for exampleFileName in os.listdir(examplePath):
      if any([exampleFileName.endswith(ext) for ext in [".cpp", ".py", ".java", ".m"]]):
        example = convertExampleSourceToDoxy(os.path.join(examplePath, exampleFileName))
        if example is not None:
          examples.append(example)

  return examples

# create a Doxyfile containing proper module paths based on Doxyfile_template
def createDoxyfile(modules):
  inputPaths = "INPUT ="
  excludePaths = "EXCLUDE ="
  examplePaths = "EXAMPLE_PATH ="
  imagePaths = "IMAGE_PATH ="

  for moduleName in modules:
    inputPath = moduleName + "/"
    examplePath = os.path.join(moduleName, "examples")
    snippetPath = os.path.join(moduleName, "doc", "doxygen", "snippets")
    testPath = os.path.join(moduleName, "tests")
    imagePath = os.path.join(moduleName, "doc", "doxygen", "images")

    if os.path.exists(os.path.join(os.getcwd(), inputPath)):
      inputPaths += " " + inputPath
    if os.path.exists(os.path.join(os.getcwd(), examplePath)):
      examplePaths += " " + examplePath + " " + snippetPath
      excludePaths += " " + examplePath + " " + snippetPath
    if os.path.exists(os.path.join(os.getcwd(), testPath)):
      excludePaths += " " + testPath
    if os.path.exists(os.path.join(os.getcwd(), imagePath)):
      imagePaths += " " + imagePath

  with open("Doxyfile_template", "r") as f:
    template = f.read()

  with open("Doxyfile", "w") as f:
    for line in template.splitlines():
      if re.match(r"INPUT  .*", line):
        f.write(inputPaths + "\n")
      elif re.match(r"EXCLUDE  .*", line):
        f.write(excludePaths + "\n")
      elif re.match(r"EXAMPLE_PATH  .*", line):
        f.write(examplePaths + "\n")
      elif re.match(r"IMAGE_PATH  .*", line):
        f.write(imagePaths + "\n")
      else:
        f.write(line + "\n")

# write language example file (e.g., examples_cpp.doxy)
def createLanguageExampleDoxy(examples):
  for language in ["cpp", "py", "java", "m"]:
    examplesInLanguage = [example for example in examples if example["language"] == language]

    # create examples menu page
    with open(os.path.join(
        "base", "doc", "doxygen", "examples_{}.doxy".format(language)),
        "w") as examplesFile:
      examplesFile.write("/**\n")
      languageName = {"cpp" : "C++", "py" : "Python", "java" : "Java", "m" : "MATLAB"}[language]
      examplesFile.write("@page examples_{} {} Examples\n".format(language, languageName))
      examplesFile.write("This is a list of all {} examples.\n".format(languageName))
      examplesFile.write("All examples can be found in the <tt>MODULE_NAME/example/</tt> " +
                         "directories.\n")
      if language == "cpp":
        examplesFile.write("Note that SCons automatically compiles (but not runs) " +
                           "all C++ examples on each run. The executables can be found in the " +
                           "same directory in which the examples reside and can be run " +
                           "directly, if <tt>LD_LIBRARY_PATH</tt> (on Linux/Mac) or " +
                           "<tt>PATH</tt> (on Windows) is set correctly.\n")

      # split examples of current language into the different modules
      moduleNames = sorted(list(set([example["moduleName"] for example in examplesInLanguage])))
      for moduleName in moduleNames:
        examplesInLanguageAndModule = [example for example in examplesInLanguage
                                       if example["moduleName"] == moduleName]
        examplesInLanguageAndModule.sort(key=lambda example: example["pageName"].lower())
        examplesFile.write("\n@section examples_{1}_module_{0} Module sgpp::{0}\n\n".format(moduleName, language))
        for example in examplesInLanguageAndModule:
          examplesFile.write("- @subpage {}\n".format(example["pageName"]))

      examplesFile.write("**/\n")

# patch Doxygen navigation tree to extend specific navtree points by default
# to ensure that users see and read those pages more likely;
# has to be called after executing Doxygen
def patchNavtree(target, source, env):
  navtreePath = os.path.join("doc", "html", "navtree.js")
  with open(navtreePath, "r") as f: navtree = f.read()

  # insert patch after
  #     $(window).load(function(){
  #       navTo(o,toroot,hashUrl(),relpath);
  #       showRoot();
  #      });
  # into the initNavTree function (in navtree.js)
  # which is called by all the *.html files
  lines = navtree.splitlines()
  patched = False

  for i, line in enumerate(lines):
    if "$(window).load(function(){" in line:
      # insert the following at line i+4
      lines[i+4:i+4] = """
  // expand some navtree points by default
  $(window).load(function() {
    try {
      // "SG++" (root node)
      expandNode(o, o.node.children[0], true, true);
      // "SG++" -> "SG++ Documentation"
      expandNode(o, o.node.children[0].children[0], true, true);
      // "SG++" -> "SG++ Documentation" -> "Manual"
      expandNode(o, o.node.children[0].children[0].children[3], true, true);
      // "SG++" -> "SG++ Documentation" -> "Manual" -> "Installation and Usage"
      expandNode(o, o.node.children[0].children[0].children[3].children[1], true, true);
      // "SG++" -> "SG++ Documentation" -> "Manual" -> "Examples"
      expandNode(o, o.node.children[0].children[0].children[3].children[3], true, true);
    } catch (err) {
      // expandNode can fail if the navtree isn't yet initialized ==> try again later
      setTimeout(arguments.callee, 0);
    }
  });
""".splitlines()
      patched = True
      break

  if not patched:
    Helper.printWarning("Doxygen navtree could not be patched.")
    return

  navtree = os.linesep.join(lines)
  with open(navtreePath, "w") as f: f.write(navtree)

# create/prepare *.doxy files for Doxygen
def prepareDoxygen(modules):
  createDoxyfile(modules)
  examples = convertExampleSourcesToDoxy(modules)
  createLanguageExampleDoxy(examples)
