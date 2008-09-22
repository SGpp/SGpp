# This file is part of pysgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.

## @package converter
# @ingroup bin
# @brief Converter lets you convert raw data formats into ARFF.
# @version $CURR$
#
# Having a typical text file containing data with one column for each
# attribute, separated by whitespaces or some other character, this tool lets
# you convert the data into the Weka-ARFF format. It allows for normalizing
# both the data dimensions to [0,1]^d and the class attribute to {-1,+1}.

from optparse import OptionParser
import sys
from tools import *

parser = OptionParser()
parser.add_option("-i", "--infile", action="append", type="string", dest="infiles", help="Specifies the inputfile to normalize. If multiple inputfiles are given they are normalized together, not separately.")
parser.add_option("-t", "--types", action="append", type="string", dest="types", help="Specifies the type of the inputfiles. Availible types are \"arrf\" and \"simple\". Can be specified for each inputfile.")
parser.add_option("-m", "--merge", action="store_true", default=False, dest="merge", help="If this is enabled, then all inputfiles are normalized and merged to one outputfile.")
parser.add_option("-n", "--nonormalization", action="store_true", default=False, dest="nonormalization", help="If this option is set, neither data (-b) nor classes (-c) are normalized, but all data is just being converted.")
parser.add_option("-b", "--border", action="store", type="float", dest="border", help="Specifies the border for the Dataset. If not set 0.05 is used.")
parser.add_option("-c", "--class", action="store", type="float", dest="c_border",metavar="BORDER", help="Specifies the classification border on which the classification data is put into different classes. If not set classvalue>=0.5 is used.")
parser.add_option("--class_min", action="store", type="float", dest="c_border_min", default=-sys.maxint-1, metavar="BORDERMIN", help="Specifies a second classification border on which the classification data is put into different classes (additionally check for classvalue<=class_min). If not set, it is ignored.")
parser.add_option("-C", "--noclasses", action="store_true", default=False, dest="noclasses", help="If this is enabled, then inputfiles have no classes.")
parser.add_option("-o", "--output", action="append", type="string", dest="outfiles", help="Specifies the output file. Can be specified for each inputfile. If not applicated, .arff is appended to each file.")
parser.add_option("--delimiter", action="store", type="string", default="", dest="delimiter", help="The delimiter separating the columns for the simple-format. Default: Split for whitespaces.")
parser.add_option("--maple", action="store_true", default=False, dest="maple", help="If enabled, write Maple-readable format.")
(options,args)=parser.parse_args()

if options.infiles == None:
	print("No inputfiles. Aborting... Help with -h")
	sys.exit(1)
	
if options.border == None:
	options.border = 0.05
	
if options.c_border == None:
	options.c_border = 0.5

if len(args) != 0:
	print("Warning: There were command-line args. Maybe forgotten to put an -i in front of filename?")

if options.outfiles != None and len(options.outfiles) > len(options.infiles):
	print("Warning: There are more outputfiles than inputfiles. Thats nonsens!?! Rest will be ignored!")
	options.outfiles = options.outfiles[:len(options.infiles)]

if options.types != None and len(options.types) > len(options.infiles):
	print("Warning: There are more types specified than inputfiles. Thats nonsens!?! Rest will be ignored!")
	options.types = options.types[:len(options.infiles)]

if options.outfiles == None:
	options.outfiles = []

if options.types == None:
	options.types = []

# Nun sind alle Arrays angelegt, und nur noch so lang wie infile[]

for i in xrange(len(options.infiles)):
	if i >= len(options.outfiles):
		options.outfiles.append(options.infiles[i] + ".arff")

	if i >= len(options.types):
		filename = options.infiles[i].lower()
		if filename.endswith("arff"):
			options.types.append("arff")
		else:
			options.types.append("simple")
	else:
		options.types[i] = options.types[i].lower()
		if not (options.types[i] == "arff" or options.types[i] == "simple"):
			print("Filetype " + options.types[i] + " is unknown. Aborting ...")
			sys.exit(1)

# Ok, alles zusammen, kann losgehen: in options.infiles[] sind alle inputfiles, 
# in options.outfiles[] sind die ausgaben, in options.merge ist ein boolean, 
# options.types[] sind alle Typen drinnen, alle Arrays sind gleich lang.

data = []

for i in xrange(len(options.infiles)):
	try:
		if options.types[i] == "arff":
			data.append(readDataARFF(options.infiles[i]))
			data[i]["filename"] = options.outfiles[i]
		elif options.types[i] == "simple":
			data.append(readDataTrivial(options.infiles[i], delim=options.delimiter))
			data[i]["filename"] = options.outfiles[i]
	except Exception:
		print("Error while reading "  + options.infiles[i] +"! Aborting...");
		sys.exit(1)
				
checkData(data)

if not options.nonormalization:
	normalize(data, options.border)
	normalizeClasses(data, options.c_border, options.c_border_min)
if not options.maple:
	writeDataARFF(data, options.merge)
else:
	writeDataMaple(data, options.merge)
