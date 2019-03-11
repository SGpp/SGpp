# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python

from builtins import str
from builtins import range
import re, sys, os

if len(sys.argv) <= 1:
    print("""Usage: python %s file [files]
Replaces all occuring old includes by the new ones in one or more files""" % (sys.argv[0]))
    print("""Recursive call, e.g. all files recursively from current directory src/sgpp via 
python ../../tmp/rename_headers.py `find -name "*" -a -not -iwholename "*.svn*" -printf "%p "`""")
    exit(1)

if not os.path.exists(sys.path[0]+'/move_headers.txt'):
    print("File move_headers.txt must be in same directory as script")
    exit(1)

# read in renames:
renames = open(sys.path[0]+'/move_headers.txt', 'r').readlines()
# skip all lines starting with #
for i in range(len(renames)-1,-1,-1):
    if re.match("#", renames[i]):
        print(renames[i])
        del renames[i]
# split them
renames = [x.split() for x in renames]

# read in files:
for f in sys.argv[1:]:
    print("Reading file", f)

    # read in file
    if not os.path.isfile(f):
        print("Error: Skipping", f, "which is not a file")
        continue
    fd = open(f, 'r')
    txt = fd.read()
    fd.close()
    txtold = str(txt)
    
    # do replacements:
    for rn in renames:
        if len(rn) == 2:
            rn = [os.path.relpath(l) for l in rn]
            # for xpp files
            txt = re.sub('"'+rn[0]+'"', '"'+rn[1]+'"', txt)
            # for .i files
            txt = re.sub('"src/sgpp/'+rn[0]+'"', '"src/sgpp/'+rn[1]+'"', txt)
            # for SConstructs
            regexp = re.compile('^'+rn[0]+'$', re.MULTILINE)
            txt = re.sub(regexp, rn[1], txt)

    # feedback
    if not txtold == txt:
        print("... modified")

        # write back file
        fd = open(f, 'w')
        fd.write(txt)
        fd.close()
