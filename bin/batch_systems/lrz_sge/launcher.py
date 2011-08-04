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

file = "dmc2008_training_4k.arff"
test = "dmc2008_testing_4k.arff"

level = 2

filename = "dmc4k10"

zeh = "identity"

additionalflags = "-a 10 -b -r 0.001 -v"

# start launcher.py by hand to get number of tasks
lrange = [0.1**(i/(40.0/6.0)) for i in range(-19,5+1)]

stats_dir = "stats/"

#################################################

stats_name = stats_dir + filename + ".stats"

if not os.environ.has_key("SGE_TASK_ID"):
    print lrange
    print "Tasks: ", len(lrange)

else:
    index = int(os.environ["SGE_TASK_ID"])-1

    command = "--mode test --zeh %s --data %s --test %s --stats %s -l %d -L %.10f --seed 6 %s" % (zeh, file, test, stats_name, level, lrange[index], additionalflags)
    commands = command.split()
   
    os.execvp("python", ["python", classifier] + commands)


