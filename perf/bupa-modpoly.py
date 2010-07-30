#!/usr/bin/python
# This file is part of sgpp, a set of tools utilizing sparse grids to solve numerical problems
# 
# Copyright (C) 2008  Joerg Blank (blankj@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os

classifier = "../bin/classifier.py"

file = " bupa-norm.arff"

level = 2

filename = "bupa"

zeh = "identity"

basetype = "modpoly"

additionalflags = "-a 40 -r 0.0001 -v -P 6"

# start launcher.py by hand to get number of tasks
lrange = [0.1**(i/(20.0/6.0)) for i in range(-1,15+1)]

stats_dir = "stats/"

#################################################

if not os.environ.has_key("SGE_TASK_ID"):
    print lrange
    print "Tasks: ", len(lrange)

else:
    index = int(os.environ["SGE_TASK_ID"])-1
    
    stats_name = stats_dir + filename + "_" + basetype + "_" + str(index) + ".stats"
    
    command = "--mode folds --zeh %s --data %s --base %s --stats %s -l %d -L %.10f --seed 6 %s" % (zeh, file, basetype, stats_name, level, lrange[index], additionalflags)
    commands = command.split()
   
    os.execvp("python", ["python", classifier] + commands)


