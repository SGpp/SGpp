# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python
# -*- coding: latin-1 -*

#############################################################################
                                    #
#############################################################################

#!/usr/bin/python
# -*- coding: latin-1 -*-
import sys

if len(sys.argv) != 2:
  print "Parameter missing. Specifiy either 1 or 2"

elif sys.argv[1] == "1":
  dat = "#!/bin/bash\n\n"
  s_tr = "./create_dataset.py -d %d -n 20000 -o chess_%02dD_tr.dat.gz\n"
  s_te = "./create_dataset.py -d %d -n 10000 -o chess_%02dD_te.dat.gz\n"
  for i in range(1,21):
    dat = dat+s_tr%(i, i)+s_te%(i,i)
  print dat

elif sys.argv[1] == "2":
  dat = "#!/bin/bash\n\n"
  s = "python ../../bin/converter.py -i chess_%02dD_tr.dat.gz -i chess_%02dD_te.dat.gz -t simple\n"
  for i in range(1,21):
    dat = dat+s%(i, i)
  print dat