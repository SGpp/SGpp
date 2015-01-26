// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#!/usr/bin/env python
import sys
import os

if __name__ == '__main__':
    filename = sys.argv[1]
    if os.path.isfile(filename):
        (filepath, ext) = os.path.splitext(filename)
        if ext in [".cpp", ".c", ".cc", ".h", ".hpp"]:
            os.system("astyle " + filename)
