import os
import sys
import random

def choose_working_directory(scratch_path):
    out = os.path.join(scratch_path, repr(random.randint(1,1000000)))
    while os.path.exists(out):
        out = os.path.join(scratch_path, repr(random.randint(1,1000000)))
    os.makedirs(out)
    print out