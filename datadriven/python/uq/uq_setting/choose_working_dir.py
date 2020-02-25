# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import os
import sys
import random

def choose_working_directory(scratch_path):
    out = os.path.join(scratch_path, repr(random.randint(1,1000000)))
    while os.path.exists(out):
        out = os.path.join(scratch_path, repr(random.randint(1,1000000)))
    os.makedirs(out)
    print( out )