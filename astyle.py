#!/usr/bin/env python
import sys
import os

if __name__ == '__main__':
    filename = sys.argv[1]
    print "argument: " + str(filename)
    if os.path.isfile(filename):
        (filepath, ext) = os.path.splitext(filename)
        if ext in [".cpp", ".c", ".cc", ".h", ".hpp"]:
            os.system("astyle " + filename)
        else:
            print "astyle: file was ignored since it has an extension unusual for c++"
