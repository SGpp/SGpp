#!/usr/bin/python
# -*- coding: latin-1 -*

#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
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