#!/usr/bin/python

import re, sys, os

if len(sys.argv) <= 1:
    print """Usage: python %s file [files]
Replaces all occuring old includes by the new ones in one or more files""" % (sys.argv[0])
    print """Recursive call, e.g. all files recursively from current directory src/sgpp via 
python ../../tmp/rename_headers.py `find -name "*" -a -not -iwholename "*.svn*" -printf "%p "`"""
    exit(1)

if not os.path.exists(sys.path[0]+'/move_headers.txt'):
    print "File move_headers.txt must be in same directory as script"
    exit(1)

# read in renames:
renames = open(sys.path[0]+'/move_headers.txt', 'r').readlines()
renames = map(lambda x: x.split(), renames)

# read in files:
for f in sys.argv[1:]:
    print "Reading file", f

    # read in file
    if not os.path.isfile(f):
        print "Error: Skipping", f, "which is not a file"
        continue
    fd = open(f, 'r')
    txt = fd.read()
    fd.close()
    txtold = str(txt)
    
    # do replacements:
    for rn in renames:
        if len(rn) == 2:
            txt = re.sub('"'+rn[0][2:]+'"', '"'+rn[1][2:]+'"', txt)
    
    # feedback
    if not txtold == txt:
        print "... modified"

        # write back file
        fd = open(f, 'w')
        fd.write(txt)
        fd.close()
