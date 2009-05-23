# This file is part of sgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## @package tools
# @ingroup bin
# @brief A collection of helper functions.
# @version $CURR$

#data = [{"data":[[x],[y]], "classes":[c]}, {"data":[[x],[y]], "classes":[c]}]
import sys, re, time, fcntl, os, random, gzip, math
from pysgpp import *

# constants
ARFF = 1
SIMPLE = 0
NOTAFILE = -1

#-------------------------------------------------------------------------------
## Opens a file. If the file ends with ".gz", automatically gzip compression
# is used for the file. Returns the filedescriptor
# @param filename file's filename
# @param mode default: "r" for read only
# @return file descriptor
def gzOpen(filename, mode="r"):
    # gzip-file
    if re.match(".*\.gz$", filename):
        # mode set for binary data?
        if not mode[-1] == "b":
            mode += "b"
        fd = gzip.open(filename, mode)
    # non gzip-file
    else:
        fd = open(filename, mode)
    return fd

#-------------------------------------------------------------------------------
## Writes a String txt to File filename, appends by default.
# Uses secure writing, i.e. locks file.
# On Windows concurrent access raises an error wich is handled.
# On Linux/Unix it should block until lock released!!
#     param: filename
#     param: txt
#     param: mode, default: "a"
def writeLockFile(filename, txt, mode="a"):
    try:
        f = gzOpen(filename, mode)
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(txt)
        fcntl.flock(f, fcntl.LOCK_UN)
        f.close()
    except IOError, e:
        sys.stderr.write("Unable to acquire lock "+str(e)+os.linesep)
        fname = filename+"."+str(time.time())
        sys.stderr.write("Writing to "+fname+os.linesep)
        f = gzOpen(fname, mode)
        f.write(txt)
        f.close()

#-------------------------------------------------------------------------------
## checks whether a file given by a filename is an ARFF-file
#     param: filename
#     returns: ARFF, SIMPLE or NOTAFILE
def isARFFFile(filename):
    try:
        # read lines until non-empty line found
        f = gzOpen(filename, 'r')
        data = f.readline().strip()
        while len(data) == 0:
            data = f.readline().strip()
        f.close()
        # check whether it starts with "@"
        if re.match('@', data):
            return ARFF
        else:
            return SIMPLE
    except:
        return NOTAFILE

#-------------------------------------------------------------------------------
## Writes String to File and checks if file existant
def writeStringToFile(s, filename, check=True):
    if check and os.path.exists(filename):
        i = raw_input("File <%s> exists. Overwrite [y/n]? " % (filename))
        if len(i) > 0 and i[0] == 'y':
            f = gzOpen(filename, 'w')
            f.write(s)
            f.close()
    else:
        f = gzOpen(filename, 'w')
        f.write(s)
        f.close()
    
#-------------------------------------------------------------------------------
## Converts a python "dataset"-structure into two object of type DataVector (X,Y)
def createDataVectorFromDataset(dataset):
    dim = len(dataset["data"])
    entries = len(dataset["data"][0])
    data = dataset["data"]
    
    x = DataVector(entries, dim)
    y = None

    for d in xrange(dim):    
        for i in xrange(entries):
            x[i*dim + d] = data[d][i]
            
    if dataset.has_key("classes"):
        classes = dataset["classes"]
        y = DataVector(entries)
        for i in xrange(entries):
            y[i] = classes[i]
            
    return (x,y)

#-------------------------------------------------------------------------------
## Converts one or two (data + optionally classes) objects of type DataVector
# to a python "dataset"-structure
def createDatasetFromDataVector(data, classes=None):
    dataset = {}
    dataset['data'] = []
    for j in range(data.getDim()):
        column = [data.get(i,j) for i in range(data.getSize())]
        dataset['data'].append(column)

    if classes:
        dataset['classes'] = [classes[i] for i in range(len(classes))]
    
    return dataset

def readDataTrivial(filename, delim = "", hasclass = True):
    fin = gzOpen(filename, "r")
    data = []
    classes = []
    for line in fin:
        sline = line.strip()
        if sline.startswith("#") or len(sline) == 0:
            continue
        
        values = []
        
        if len(delim) == 0:
            values = sline.split()
        else:
            values = sline.split(delim)
            values = map(lambda x: x.strip(), values)
            values = filter(lambda x: len(x) > 0, values)
            
        if len(data) == 0:
            if hasclass:
                data = [[] for i in range(len(values) - 1)] 
            else:
                data = [[] for i in range(len(values))] 
        
        if hasclass:
            for i in xrange(len(values)-1):
                data[i].append(float(values[i]))
                
            classes.append(float(values[-1]))
        else:
            for i in xrange(len(values)):
                data[i].append(float(values[i]))


    fin.close()
    
    if hasclass:
        return {"data":data, "classes":classes, "filename":filename}
    else:
        return {"data":data, "filename":filename}

## Reads in an ARFF file:
# The data is stored in lists. There is a value list for every dimension of the data set. e.g. 
# [[2, 3],[1, 1]] are the data points P_1(2,1) and P_2(3,1)
#
# @param filename the file's filename that should be read
# @return returns a set of a array with the data (named data), a array with the classes (named classes) and the filename named as filename
def readDataARFF(filename):
    fin = gzOpen(filename, "r")
    data = []
    classes = []
    hasclass = False

    # get the different section of ARFF-File
    for line in fin:
        sline = line.strip().lower()
        if sline.startswith("%") or len(sline) == 0:
            continue

        if sline.startswith("@data"):
            break
        
        if sline.startswith("@attribute"):
            value = sline.split()
            if value[1].startswith("class"):
                hasclass = True
            else:
                data.append([])
    
    #read in the data stored in the ARFF file
    for line in fin:
        sline = line.strip()
        if sline.startswith("%") or len(sline) == 0:
            continue

        values = sline.split(",")
        if hasclass:
            classes.append(float(values[-1]))
            values = values[:-1]
        for i in xrange(len(values)):
            data[i].append(float(values[i]))
            
    # cleaning up and return
    fin.close()
    return {"data":data, "classes":classes, "filename":filename}


#-------------------------------------------------------------------------------
## Writes gnuplot data of function into file
# 
def writeGnuplot(filename, grid, alpha, resolution):
    p = DataVector(1,2)
    fout = file(filename, "w")

    for x in xrange(resolution):
        for y in xrange(resolution):
            p[0] = float(x) / (resolution - 1)
            p[1] = float(y) / (resolution - 1)
            pc = grid.createOperationEval().eval(alpha, p)
            fout.write("%f %f %f\n" % (p[0], p[1], pc))
        fout.write("\n")
    fout.close()
    return

#-------------------------------------------------------------------------------
## Writes gnuplot data of grid into file
# 
def writeGnuplotGrid(filename, grid):
    p = DataVector(1,2)
    fout = file(filename, "w")
    s = grid.__str__()
    s = s.split("], [")
    for gp in s:
        (l1,i1,l2,i2) = re.search("(\d+)[, ]*(\d+)[, ]*(\d+)[, ]*(\d+)", gp).groups()
        fout.write("%f %f\n" % (int(i1)*2**(-int(l1)), int(i2)*2**(-int(l2))))
    fout.close()
    return

def writeDataARFF(data, merge=False):
    if len(data) == 0:
        return
    # is data a dataset or is it a list of datasets?
    if isinstance(data, dict):
        data = [data]
    
    hasHeader = False
    for dataset in data:
    	if hasHeader == False or merge == False:
            hasHeader = True
            fout = gzOpen(dataset["filename"], "w")
            fout.write("@RELATION \"%s\"\n\n" % dataset["filename"])
            fstring = ""

            if isinstance(dataset["data"], DataVector):
                dim = dataset["data"].getDim()
            else:
                dim = len(dataset["data"])
                   
            for i in xrange(dim):
                fout.write("@ATTRIBUTE x%d NUMERIC\n" % i)
                fstring = fstring + "%s,"
           
            hasclass = False
            if dataset.has_key("classes"):
           	hasclass = True
           	fout.write("@ATTRIBUTE class NUMERIC\n")
           	fstring = fstring + "%s"
            else:
           	fstring = fstring.strip(',')
                
            fstring = fstring + "\n"
            fout.write("\n@DATA\n")
        
        if isinstance(dataset["data"], DataVector):
            num_rows = dataset["data"].getSize()
            for row in xrange(num_rows):
                lout = []
                for column in xrange(dim):
                    lout.append(dataset["data"][row*dim+column])
                if hasclass:
                    lout.append(dataset["classes"][row])
                fout.write(fstring % tuple(lout))
        else:
            num_rows = len(dataset["data"][0])
            for row in xrange(num_rows):
                lout = []
                for column in xrange(dim):
                    lout.append(dataset["data"][column][row])
                if hasclass:
                    lout.append(dataset["classes"][row])
                fout.write(fstring % tuple(lout))
        
        if merge == False:
        	fout.close()
    
    if merge == True:
        fout.close()
    return


#-------------------------------------------------------------------------------
## Writes a data object to a file, specified by data["filename"]+".maple"
# @param data a data object
# @param merge set to True, iff data is a list of data objects and they
#        should be joined
def writeDataMaple(data, merge):
    if merge == True:
        print "Ignoring merge. Not implemented yet!"
    
    if len(data) == 0:
        return
    
    for dataset in data:
        s = ""

        # transform
        lclass = []
        for row in xrange(len(dataset["data"][0])):
            lout = []
            for column in xrange(len(dataset["data"])):
                lout.append(dataset["data"][column][row])
            if dataset.has_key("classes"):
                lclass.append(dataset["classes"][row])
            s += '['+','.join(map(lambda l: str(l), lout))+'],\n'
        s = "X := Matrix([%s]);\n" %(s[:-2])
        s += "Y := Vector([%s]);\n" % (','.join(map(lambda l: str(l), lclass)))

        # write out 
        fout = open(dataset["filename"]+".maple", "w")
        fout.write(s)
        fout.close()

    return

#-------------------------------------------------------------------------------
## Writes a DataVector object to a file, specified by filename.
# Output in the file is "X := Matrix([[...],[...],...,[...]]);"
# @param data a DataVector object
# @param filename the file's name
# @param format (optional) format specifier, default: "%s"
# @param maple_name (optional) name of variable in Maple, default: "X"
# @param check (optional) if set to true, program will ask before overwriting files
def writeDataVectorMaple(data, filename, format="%s", maple_name="X", check=True):
    numrows = data.getSize()
    numcols = data.getDim()
    s = "%s := Matrix([\n" % (maple_name)
    for row in xrange(numrows):
        col_list = []
        for col in xrange(numcols):
            col_list.append(format % (data.get(row, col)))
        s += "["+",".join(col_list)+"]"
        if row < numrows-1:
            s += ","
        s += "\n"
    s += "], datatype=float);"
    writeStringToFile(s, filename, check=check)


#-------------------------------------------------------------------------------
## Writes information that is needed for the normalization of data to a file.
# Using this information one can then later on reverse the normalization or
# normalize further data.
# @param filename a filename
# @param border offset for normalization
# @param minvals the (original) minimum value of each attribute
# @param maxvals the (original) maximum value of each attribute
def writeNormfile(filename, border, minvals, maxvals):
    s  = "border: %s\n" % (border)
    s += "min:    %s\n" % (' '.join(map(lambda l: "%s"%l, minvals)))
    s += "max:    %s\n" % (' '.join(map(lambda l: "%s"%l, maxvals)))
    writeStringToFile(s, filename)
    
#-------------------------------------------------------------------------------
## Reads information that is needed for the normalization of data from a file.
# @param filename a filename
# @return (border, minvals, maxvals, deltavals)
# border: offset for normalization
# minvals: the minimum value of each attribute
# maxvals: the maximum value of each attribute
# deltavals: (max-min)/(1.0-2*border), provided for convenience
def readNormfile(filename):
    fd = gzOpen(filename, 'r')
    data = fd.readlines()
    fd.close()
    try:
        border = float((data[0].strip()).split(None)[1])
        minvals = map(lambda l: float(l), (data[1].strip()).split(None)[1:])
        maxvals = map(lambda l: float(l), (data[2].strip()).split(None)[1:])
        deltavals = map(lambda x,y: ((y-x)/(1.0-2.0*border)), minvals, maxvals)
    except:
        raise Exception("ERROR: Unable to read \"%s\"\n" % (filename))
    return (border, minvals, maxvals, deltavals)


#-------------------------------------------------------------------------------
## Normalize values of input vectors on the segment [0,1]
# @param data Dataset
# @param border Specifies border of the dataset, will be added to the normalized value 
# @param filename Filename of normfile (optional)
# @param minvals Array of normalization boundary min values (one per dimension) (optional)
# @param maxvals Array of normalization boundary max values (one per dimension) (optional)
def normalize(data, border=0.0, filename=None, minvals=None, maxvals=None):
    # check parameters
    if len(data) == 0:
        raise ValueError, "Wrong or no data."
    if minvals and maxvals:
        if (len(minvals) <> len(data[0]["data"]) or 
            len(maxvals) <> len(data[0]["data"])):
            raise ValueError, "Wrong number of min- or max-values."
        lmin = minvals
        lmax = maxvals
    else:
        lmin = []
        lmax = []

        for datadim in data[0]["data"]:
            lmin.append(datadim[0])
            lmax.append(datadim[0])

        for dataset in data:
            for dim in xrange(len(dataset["data"])):
                cmin = min(dataset["data"][dim])
                lmin[dim] = min(cmin, lmin[dim])

                cmax = max(dataset["data"][dim])
                lmax[dim] = max(cmax, lmax[dim])

    # delta values
    ldelta = map(lambda x,y: ((y-x)/(1.0-2.0*border)), lmin, lmax)

    # write normalization data to file:
    if filename:
        writeNormfile(filename, border, lmin, lmax, ldelta)
    
    for dataset in data:
        for dim in xrange(len(dataset["data"])):
            # special handling for the case that all max==min, i.e. all
            # attribute values are equal: set to 0.5
            if ldelta[dim] == 0:
                dataset["data"][dim] = map(lambda x: 0.5,dataset["data"][dim])
            else:
                dataset["data"][dim] = map(lambda x: (x-lmin[dim]) / ldelta[dim] + border,dataset["data"][dim])
    return


#-------------------------------------------------------------------------------
## Divides the class values in two categories
# @param data Dataset 
# @param border Classes will be differentiated between greater and less then border 
# @param minborder All classes under the minborder are processes as they are over border
def normalizeClasses(data, border=0.0, minborder=-sys.maxint-1):
    def separate(x):
        if x >= border or x < minborder:
            return 1
        else:
            return -1
        
    for dataset in data:
        dataset["classes"] = map(separate, dataset["classes"])
    return

#-------------------------------------------------------------------------------
## Validates Dataset
# @param data Dataset 
def checkData(data):
    if len(data) == 0:
    	print("No data loaded. Aborting!")
    	sys.exit(1)

    hasClasses = False
    
    data_len = len(data[0]["data"])
    if data[0].has_key("classes"):
    	hasClasses = True
    
    for dataset in data:
    	if data_len != len(dataset["data"]):
    		print("Error! Can't merge data due to different amount of dimensions!")
    		sys.exit(1)
    	if dataset.has_key("classes") != hasClasses:
    		print("Error! Can't merge data, because some files have classification data and some not!")
    		sys.exit(1)
    		
    return


#-------------------------------------------------------------------------------
# perform split of dataset in num_partitions partitions for n-fold-cv
# split data into folds,
# return ([data1,data2,...,datan], [classes1,classes2,...,classesn])
#-------------------------------------------------------------------------------
def split_n_folds(data, num_partitions, seed=None):
    dim = len(data["data"])
    size = len(data["data"][0])

    random.seed(seed)
    seq = range(size)
    random.shuffle(seq)

    dvec = []
    cvec = []

    size_left = size
    index = 0

    for i in xrange(num_partitions):
        size_fold = size_left/(num_partitions-i)
        dvec.append(DataVector(size_fold, dim))
        cvec.append(DataVector(size_fold))
        for element in xrange(size_fold):
            for d in xrange(dim):
                dvec[i][element*dim + d] = data["data"][d][seq[index]]
            #@todo: this doesn't work for regression, because the last parameter is not necessary class
            cvec[i][element] = data["classes"][seq[index]]
            index += 1
        size_left = size_left-size_fold
    
    return (dvec, cvec)

#-------------------------------------------------------------------------------
# perform sequential(!) split of dataset in num_partitions partitions for n-fold-cv
# split data into folds,
# return ([data1,data2,...,datan], [classes1,classes2,...,classesn])
#-------------------------------------------------------------------------------
def split_n_folds_sequential(data, num_partitions):
    dim = len(data["data"])
    size = len(data["data"][0])

    dvec = []
    cvec = []
    size_left = size
    index = 0

    for i in xrange(num_partitions):
        size_fold = size_left/(num_partitions-i)
        dvec.append(DataVector(size_fold, dim))
        cvec.append(DataVector(size_fold))
        for element in xrange(size_fold):
            for d in xrange(dim):
                dvec[i][element*dim + d] = data["data"][d][index]
            cvec[i][element] = data["classes"][index]
            index += 1
        size_left = size_left-size_fold

    return (dvec, cvec)


#-------------------------------------------------------------------------------
# perform split of dataset in num_partitions partitions for stratified n-fold-cv
# split data into folds,
# return ([data1,data2,...,datan], [classes1,classes2,...,classesn])
#-------------------------------------------------------------------------------
def split_n_folds_stratified(data, num_partitions, seed=None):
    dim = len(data["data"])
    size = len(data["data"][0])
    # split in pos. and neg sets
    neg = []
    pos = []
    for i in xrange(size):
        if data["classes"][i] < 0:
            neg.append(i)
        else:
            pos.append(i)
    
    # shuffel using seed
    random.seed(seed)
    random.shuffle(neg)
    random.shuffle(pos)

    dvec = []
    cvec = []

    size_left_pos = len(pos)
    size_left_neg = len(neg)
    index_pos = 0
    index_neg = 0

    for i in xrange(num_partitions):
        size_fold_pos = size_left_pos/(num_partitions-i)
        size_fold_neg = size_left_neg/(num_partitions-i)
        dvec.append(DataVector(size_fold_pos+size_fold_neg, dim))
        cvec.append(DataVector(size_fold_pos+size_fold_neg))
        # add data with class pos first
        for element in xrange(size_fold_pos):
            for d in xrange(dim):
                dvec[i][element*dim + d] = data["data"][d][pos[index_pos]]
            cvec[i][element] = data["classes"][pos[index_pos]]
            index_pos += 1
        size_left_pos = size_left_pos-size_fold_pos
        # then add data with class neg
        for element in xrange(size_fold_neg):
            for d in xrange(dim):
                dvec[i][(size_fold_pos+element)*dim + d] = data["data"][d][neg[index_neg]]
            cvec[i][size_fold_pos+element] = data["classes"][neg[index_neg]]
            index_neg += 1
        size_left_neg = size_left_neg-size_fold_neg

    return (dvec, cvec)

#-------------------------------------------------------------------------------
## perform sequential split of a DataVector into two DataVectors
# @param data DataVector to split
# @param proportion split into proportion, (1-proportion)
# @return (DataVector1, DataVector2)
#-------------------------------------------------------------------------------
def split_DataVector_by_proportion(data, proportion):
    dim = data.getDim()
    size = data.getSize()

    splitpoint = int(math.floor(size*proportion))
    dv1 = DataVector(splitpoint, dim)
    dv2 = DataVector(size-splitpoint, dim)
    row = DataVector(1,dim)
    
    # copy
    for i in xrange(splitpoint):
        data.getRow(i, row)
        dv1.setRow(i, row)
    for i in xrange(size-splitpoint):
        data.getRow(i+splitpoint, row)
        dv2.setRow(i, row)

    if len(dv1)+len(dv2) <> len(data): raise Exception("data length doesn't match")
    return (dv1, dv2)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## An array containing all modes and descriptions
zeh_modes = {
    "laplace" : "Classical Laplacian. See OpLaplaceAdaptive",
    "identity" : "Identity matrix, most efficient. See OpPseudo",
    "ratio" : "Preferres quadratical supports. See OpPseudo",
    "levelsum" : "Sum of the levels, scaled by the gridlevel (usually 2 for adaptive SGs). See OpPseudo",
    "energy" : "Energy-norm-like SGs. See OpPseudo",
    "copy" : "Sum of the rows of classical Laplacian. See OpLaplaceAdaptive",
    "pseudounit" : "???"}
    
## base function types
base_types = {
#    "linear" : {"base" : SLinearBase, "b" : SGridOperationB, "t" : test_dataset_linear, "laplace" : SGridOperationLaplace},
#    "modlinear" : {"base" : SModLinearBase, "b" : SGridModOperationB, "t" : test_dataset_modlin},
#    "poly" : {"base" : "SPolyBase", },
              }


## Class Matrix that incorporates settings and actions for applying
# the matrix C and computing the RHS b.
# @todo fully update to pysgpp
class Matrix:
    def __init__(self, grid, x, l, mode, base = None):
        self.grid = grid
        self.x = x
        self.l = l
        self.B = grid.createOperationB()
        #val: base is obviously not used
        #self.base = base
        self.CMode = mode.lower()
        
        if self.CMode == "laplace":
            self.C = grid.createOperationLaplace()
        elif self.CMode == "identity":
            pass
        elif self.CMode == "ratio":
            self.C = OpPseudo(grid)
            self.C.initRatio()
        elif self.CMode == "levelsum":
            self.C = OpPseudo(grid)
            self.C.initLevelSum()
        elif self.CMode == "energy":
            self.C = OpPseudo(grid)
            self.C.initEnergy()
        elif self.CMode == "copy":
            self.C = OpPseudo(grid)
            self.C.initCopyFrom()
        elif self.CMode == "pseudounit":
            self.C = OpPseudo(grid)
            self.C.initPseudoUnit()
        
    
    def generateb(self, y):
        b = DataVector(self.grid.getStorage().size())
        self.B.mult(y, self.x, b)
        return b
    
    def ApplyMatrix(self, alpha, result):
        temp = DataVector(self.x.getSize())
        M = self.x.getSize();
    
        self.B.multTranspose(alpha, self.x, temp)
        self.B.mult(temp, self.x, result)
        
        temp = DataVector(alpha.getSize())

        if self.CMode == "laplace":
            # @todo: implement
            self.C.mult(alpha, temp)
            result.axpy(M*self.l, temp)

        elif self.CMode == "identity":
            result.axpy(M*self.l, alpha)
            
        elif self.CMode == "ratio":
            # @todo: implement
            self.C.applyRatio(alpha, temp)
            result.axpy(M*self.l, temp)
            
        elif self.CMode == "levelsum":
            # @todo: implement
            self.C.applyRatio(alpha, temp)
            result.axpy(M*self.l, temp)

        elif self.CMode == "energy":
            # @todo: implement
            self.C.applyRatio(alpha, temp)
            result.axpy(M*self.l, temp)

        elif self.CMode == "copy":
            # @todo: implement
            self.C.applyRatio(alpha, temp)
            result.axpy(M*self.l, temp)
        
        elif self.CMode == "pseudounit":
            # @todo: implement
            self.C.applyRatio(alpha,temp)
            result.axpy(M*self.l,temp)

#-------------------------------------------------------------------------------
# saves/restores grid from file
# set appropriate modes for automatic destructions
#-------------------------------------------------------------------------------

def restoreGrid(text):
    #@todo: is there any control of correctness of text?
    return Grid.unserialize(text)


def saveGrid(grid):
    return grid.serialize()

#-------------------------------------------------------------------------------
# read/write alpha
#-------------------------------------------------------------------------------


def writeAlphaARFF(filename, alpha):
    fout = gzOpen(filename, "w")
    fout.write("@RELATION \"%s ALPHAFILE\"\n\n" % filename)
    fout.write("@ATTRIBUTE alpha NUMERIC\n")
    
    fout.write("\n@DATA\n")
    
    for i in xrange(len(alpha)):
        fout.write("%1.20f\n" % alpha[i])
    
    fout.close()

def readAlphaARFF(filename):
    try:
        data = readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        sys.exit(1)
    
    alpha = DataVector(len(data["data"][0]), 1)
    
    for i in xrange(len(data["data"][0])):
        alpha[i] = data["data"][0][i]
    
    return alpha

def writeGrid(filename, grid):
    text = saveGrid(grid)
    fout = gzOpen(filename, "w")
    fout.write(text)
    fout.close()

def readGrid(filename):
    fin = gzOpen(filename, "r")
    text = fin.read()
    fin.close()
    
    return restoreGrid(text)

def writeCheckpoint(filename, grid, alpha, adaption = None, fold = None):
    adapt_str = ""
    fold_str = ""
    if adaption != None:
        adapt_str = ".a%d" % (adaption)
    if fold != None:
        fold_str = ".f%d" % (fold)
    writeAlphaARFF("%s%s%s.alpha.arff.gz" % (filename, fold_str, adapt_str), alpha)
    writeGrid("%s%s%s.grid.gz" % (filename, fold_str, adapt_str), grid)
    
def writeStats(filename, txt, mode = "a"):
    writeLockFile(filename + ".stats", txt, mode)

def readCheckpoint(filename):
    alpha = readAlphaARFF(filename+".alpha.arff")
    grid = readGrid(filename+".grid")
    
    return grid, alpha

