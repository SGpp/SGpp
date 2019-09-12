#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import re
import sys



def processFile(path, performFix):
  path = os.path.abspath(path)
  if os.path.splitext(path)[1] not in [".cpp", ".hpp"]: return True
  with open(path, "r") as f: source = f.read()

  # when checking examples or tests, there might be includes of helper
  # headers that are in the examples or tests directory
  # => in this case, only search for includes that contain a sgpp directory
  regex = (r"^#include \"(.*)\"$" if (os.sep + "src" + os.sep) in path else
           r"^#include \"(.*/sgpp/.*)\"$")
  matches = list(re.finditer(regex, source, flags=re.MULTILINE))
  if len(matches) == 0: return True

  if performFix:
    print("{}:".format(path))

    # for all relative include lines in reverse order,
    # since we are modifying the source on-the-fly based on indices
    for match in matches[::-1]:
      includePath = match.group(1)

      # if the include path starts with sgpp/ or if we include a library
      # like omp.h, cl.h, etc., # then just change the "" to <>
      if ((os.path.splitext(includePath)[1] == ".hpp") and
          (not includePath.startswith("sgpp/"))):
        includePath = convertIncludePath(path, includePath)

      oldIncludeLine = match.group()
      newIncludeLine = "#include <{}>".format(includePath)
      print(oldIncludeLine)
      print("-> {}".format(newIncludeLine))

      source = (source[:match.start()] +
          newIncludeLine + source[match.end():])

    print("")
    with open(path, "w") as f: f.write(source)
  else:
    print(("{}:0: warning: #include with \"quotation marks\" found. "
           "Use <angle brackets> instead and do not use relative "
           "paths.  [build/include] [5]").format(path), file=sys.stderr)
    return False



def convertIncludePath(sourcePath, includePath):
  includePath = os.path.join(os.path.dirname(sourcePath), includePath)
  includePath = os.path.normpath(includePath)
  components = includePath.split(os.sep)
  assert (components[2] == "sgpp") or components[2].startswith("sgpp_")
  includePath = os.path.join(*components[2:])
  return includePath



def main():
  parser = argparse.ArgumentParser(description="Checks for relative #include.")
  parser.add_argument("-r", "--recursive", action="store_true",
                      help="Recursively check a directory instead of a single file.")
  parser.add_argument("-f", "--fix", action="store_true",
                      help="Try to fix the relative #includes in-place, "
                           "by converting them to absolute #includes.")
  parser.add_argument("path", metavar="PATH",
                      help="Path to file or directory (with -r)")
  args = parser.parse_args()

  if args.recursive:
    success = True

    # for all source/header files
    for root, dirs, files in os.walk(args.path):
      for file_ in files:
        success = processFile(os.path.join(root, file_), args.fix) and success
  else:
    success = processFile(args.path, args.fix)

  sys.exit(0 if success else 1)



if __name__ == "__main__":
  main()
