#!/usr/bin/python
# -*- coding: latin-1 -*-

dat = "#!/bin/bash\n\n"
#s_tr = "./create_dataset.py -d %d -n 1000 -o cube_%02dD_tr.dat\n"
#s_te = "./create_dataset.py -d %d -n 500 -o cube_%02dD_te.dat\n"
#for i in range(60):
#  dat = dat+s_tr%(i, i)+s_te%(i,i)
s = "python ~/workspace/SGClass/trunk/dataset/converter.py -i cube_%02dD_tr.dat -i cube_%02dD_te.dat -t simple -b 0.015625\n"
for i in range(1,61):
  dat = dat+s%(i, i)

print dat
