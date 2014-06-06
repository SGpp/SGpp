#!/usr/bin/python
# -*- coding: latin-1 -*-
import sys

if len(sys.argv) != 2:
  print "Parameter missing. Specifiy either 1 or 2"

elif sys.argv[1] == "1":
  dat = "#!/bin/bash\n\n"
  s_tr = "./create_dataset.py -d %d -n 1000 -o cube_%02dD_tr.dat.gz\n"
  s_te = "./create_dataset.py -d %d -n 500 -o cube_%02dD_te.dat.gz\n"
  for i in range(1,61):
    dat = dat+s_tr%(i, i)+s_te%(i,i)
  print dat

elif sys.argv[1] == "2":
  dat = "#!/bin/bash\n\n"
  s = "python ../../bin/converter.py -i cube_%02dD_tr.dat.gz -i cube_%02dD_te.dat.gz -t simple -b 0.015625\n"
  for i in range(1,61):
    dat = dat+s%(i, i)
  print dat
