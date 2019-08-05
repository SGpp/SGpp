#!/usr/bin/env python

import argparse
import os
import re



def convertIncludePath(sourcePath, includePath):
  includePath = os.path.join(os.path.dirname(sourcePath), includePath)
  includePath = os.path.normpath(includePath)
  components = includePath.split(os.sep)
  assert (components[2] == "sgpp") or components[2].startswith("sgpp_")
  includePath = os.path.join(*components[2:])
  return includePath



def main():
  parser = argparse.ArgumentParser(description=
      "Recursively converts relative #include paths to absolute paths.")
  parser.add_argument("--directory", default=".",
                      help="Root directory of SG++ (default: \".\").")
  parser.add_argument("-n", "--dry-run", action="store_true",
                      help="Don't change anything, only show "
                           "what would change.")
  args = parser.parse_args()
  rootDir = args.directory
  
  # for all subdirectories
  for module in sorted(os.listdir(rootDir)):
    dirExists = (lambda *args:
        os.path.isdir(os.path.join(rootDir, module, *args)))
    
    # check if subdirectory corresponds to a module
    if dirExists() and dirExists("src") and dirExists("src", "sgpp"):
      moduleSourcePath = os.path.join(rootDir, module, "src")
      
      # for all source/header files
      for root, dirs, files in os.walk(moduleSourcePath):
        for file_ in files:
          if os.path.splitext(file_)[1] not in [".cpp", ".hpp"]: continue
          
          path = os.path.join(root, file_)
          with open(path, "r") as f: source = f.read()
          matches = list(re.finditer(
              r"^#include \"(.*)\"$", source, flags=re.MULTILINE))
          if len(matches) == 0: continue
          
          print("{}:".format(path))
          
          # for all relative include lines in reverse order,
          # since we are modifying the source on-the-fly based on indices
          for match in matches[::-1]:
            includePath = match.group(1)
            
            # if the include path starts with sgpp/,
            # then just change the "" to <>
            if not includePath.startswith("sgpp/"):
              includePath = convertIncludePath(path, includePath)
            
            oldIncludeLine = match.group()
            newIncludeLine = "#include <{}>".format(includePath)
            print(oldIncludeLine)
            print("-> {}".format(newIncludeLine))
            
            source = (source[:match.start()] +
                newIncludeLine + source[match.end():])
          
          print("")
          
          if not args.dry_run:
            with open(path, "w") as f: f.write(source)



if __name__ == "__main__":
  main()
