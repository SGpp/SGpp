#!/usr/bin/env python
import sys
import os

if __name__ == '__main__':
    filename = sys.argv[1]
    if os.path.isfile(filename):
        (filepath, ext) = os.path.splitext(filename)
        if ext in [".cpp", ".c", ".cc", ".h", ".hpp"]:
            os.system("astyle " + filename)
