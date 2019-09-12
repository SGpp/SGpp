#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import sys



def processFile(path):
  path = os.path.abspath(path)
  if os.path.splitext(path)[1] not in [".cpp", ".hpp"]: return True
  with open(path, "r") as f: source = f.read()

  COPYRIGHT_BANNER = r"""// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

"""

  if not source.startswith(COPYRIGHT_BANNER):
    print(("{}:0: warning: No SG++ copyright message found or existing "
           "copyright message is not the standard SG++ copyright "
           "message.  [legal/copyright] [5]").format(path), file=sys.stderr)
    return False

  if any([(("Author" in x) or ("Created" in x) or
           ("author" in x) or ("created" in x))
          for x in source.lower().splitlines()[:10]]):
    print(("{}:0: warning: Author or creation date in copyright message "
           "found.  [legal/copyright] [5]").format(path), file=sys.stderr)
    return False

  return True



def main():
  parser = argparse.ArgumentParser(description="Checks copyright banners.")
  parser.add_argument("-r", "--recursive", action="store_true",
                      help="Act on a directory instead of a single file.")
  parser.add_argument("path", metavar="PATH",
                      help="Path to file or directory (with -r)")
  args = parser.parse_args()

  if args.recursive:
    success = True

    # for all source/header files
    for root, dirs, files in os.walk(args.path):
      for file_ in sorted(files):
        if os.path.splitext(file_)[1] in [".cpp", ".hpp"]:
          success = processFile(os.path.join(root, file_)) and success
  else:
    success = processFile(args.path)

  sys.exit(0 if success else 1)



if __name__ == "__main__":
  main()
