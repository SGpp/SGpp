# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##
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
## @brief A value pair is added to a dictionary's value entry. 
#
# Each entry of the 
# dictionary is a list. If the dictionary has no entry for key, [value] is 
# added. Otherwise value is appended to the list under key.
# @param dict the dictionary
# @param key the key
# @param val the value
def appendToDict(dict, key, val):
    if dict.has_key(key):
        dict[key].append(val)
    else:
        dict[key] = [val]

#-------------------------------------------------------------------------------
## @brief Opens a file. If the file ends with ".gz", automatically gzip compression
# is used for the file. 
#
# Returns the filedescriptor
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
## @brief Writes a String txt to File filename, appends by default.
# Uses secure writing, i.e. locks file.
#
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
## @brief Checks whether a file given by a filename is an ARFF-file
#  @param filename the file's name
#  @return ARFF, SIMPLE or NOTAFILE
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
## @brief Writes String to File and checks if file existant
# @param s some text
# @param filename filename (including relative or absolute path)
# @param check (optional) set to False to overwrite without checking for existance
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
## @brief Reads in (multidimensional) data from a delimiter separated data file. 
#
# Last column is assumend to contain class values, if <tt>hasclass=True</tt>.
# The data is stored in a dictionary, thus either
# {"data": DataMatrix, "classes": DataVector, "filename": filename} or
# {"data": DataMatrix, "filename": filename}
#
# @param filename the file's filename that should be read
# @param delim (optional) separator between columns. Default: whitespaces
# @param hasclass (optional) sets, whether last column contains class attribute. Default: True
# @return returns the data
def readDataTrivial(filename, delim = None, hasclass = True):
    fin = gzOpen(filename, "r")
    data = []
    classes = []
    for line in fin:
        sline = line.strip()
        # skip empty lines and comments
        if sline.startswith("#") or len(sline) == 0:
            continue

        # split and convert 
        values = sline.split(delim)
        values = map(lambda x: x.strip(), values)
        values = filter(lambda x: len(x) > 0, values)
        values = map(lambda x: float(x), values)
        
        if hasclass:
            data.append(values[:-1])
            classes.append(values[-1])
        else:
            data.append(values)

    # cleaning up and return
    fin.close()
    if hasclass:
        return {"data": DataMatrix(data), "classes": DataVector(classes), "filename":filename}
    else:
        return {"data": DataMatrix(data), "filename":filename}

#-------------------------------------------------------------------------------
## @brief Reads in (multidimensional) data from an ARFF file.
#
# The data is stored in a dictionary, thus either
# {"data": DataMatrix, "classes": DataVector, "filename": filename} or
# {"data": DataMatrix, "filename": filename}, depending whether one of the attributes is called "class[es]"
#
# @param filename the file's filename that should be red
# @return returns the data
def readDataARFF(filename):
    fin = gzOpen(filename, "r")
    data = []
    classes = []
    hasclass = False

    # get the different section of ARFF-File
    for line in fin:
        sline = line.strip().lower()
        # skip comments and empty lines
        if sline.startswith("%") or len(sline) == 0:
            continue

        if sline.startswith("@data"):
            break
        
        if sline.startswith("@attribute"):
            value = sline.split()
            if value[1].startswith("class"):
                hasclass = True
    
    #read in the data stored in the ARFF file
    for line in fin:
        sline = line.strip()
        # skip comments and empty lines
        if sline.startswith("%") or len(sline) == 0:
            continue

        # split and convert 
        values = sline.split(",")
        values = map(lambda x: x.strip(), values)
        values = filter(lambda x: len(x) > 0, values)
        values = map(lambda x: float(x), values)
        
        if hasclass:
            data.append(values[:-1])
            classes.append(values[-1])
        else:
            data.append(values)
            
    # cleaning up and return
    fin.close()
    if hasclass:
        return {"data": DataMatrix(data), "classes": DataVector(classes), "filename":filename}
    else:
        return {"data": DataMatrix(data), "filename":filename}


#-------------------------------------------------------------------------------
## @brief Opens and read the (multidimensional) data of an ARFF (or plain whitespace-separated data) file.
# Assumes that class information is available. Format is
# {"data": DataMatrix, "classes": DataVector, "filename": filename} or
#
# @param filename filename of the file
# @return the data stored in the file as a set of arrays
def readData(filename):
    try:
        if isARFFFile(filename):
            data = readDataARFF(filename)
        else:
            data = readDataTrivial(filename)
    except Exception, e:
        print ("An error occured while reading " + filename + "!")
        raise e
        
    if data.has_key("classes") == False:
        print ("No classes found in the given File " + filename + "!")
        sys.exit(1)
        
    return data

#-------------------------------------------------------------------------------
## @brief Evaluates function on a full grid in the domain, and writes evaluation points
# to a file.
# The output is suitable for Gnuplot.
#
# @param filename Filename to which data is written
# @param grid Grid
# @param alpha Corresponding coefficient DataVector
# @param resolution Number of sampling points per dimension
# @param mode {'w'|'a'} to write or append, default 'w' (optional)
# @param data points to plot (optional)
# @param fvals corresponding function values (optional)
def writeGnuplot(filename, grid, alpha, resolution, mode="w", data=None, fvals=None):
    p = DataVector(grid.getDimension())
    fout = gzOpen(filename, mode)

    # evaluate 1d function
    if grid.getDimension() == 1:
        fout.write("#set term png truecolor enhanced\n")
        fout.write("#set out '%s.png'\n" % (filename))
        if data and fvals:
            fout.write("plot '-' w p lw 2, '-' w l\n")
            for i in xrange(len(fvals)):
                fout.write("%g %g\n" % (data.get(i,0), fvals[i]))
            fout.write("e\n")
        else:
            fout.write("plot '-' w l\n")
        for x in xrange(resolution):
                p[0] = float(x) / (resolution - 1)
                pc = createOperationEval(grid).eval(alpha, p)
                fout.write("%f %f\n" % (p[0], pc))
    # evaluate 2d function
    elif grid.getDimension() == 2:
        fout.write("#set term png truecolor enhanced\n")
        fout.write("#set out '%s.png'\n" % (filename))
        if data and fvals:
            fout.write("splot '-' w p lw 2, '-' w pm3d\n")
            for i in xrange(len(fvals)):
                fout.write("%g %g %g\n" % (data.get(i,0), data.get(i, 1), fvals[i]))
            fout.write("e\n")
        else:
            fout.write("splot '-' w pm3d\n")
        for x in xrange(resolution):
            for y in xrange(resolution):
                p[0] = float(x) / (resolution - 1)
                p[1] = float(y) / (resolution - 1)
                pc = createOperationEval(grid).eval(alpha, p)
                fout.write("%f %f %f\n" % (p[0], p[1], pc))
            fout.write("\n")
        fout.write("e\n")
    # can't plot anything else
    else:
        sys.stderr.write("Error! Can't plot grid with dimensionality %d..." % (grid.getDimension()))
    fout.close()
    return

#-------------------------------------------------------------------------------
## @brief Evaluates function on a full grid in the domain, and writes evaluation points
# to a file.
# The output is suitable for Gnuplot.
#
# @param filename Filename to which data is written
# @param dim dimension
# @param fctn function
# @param resolution Number of sampling points per dimension
# @param mode {'w'|'a'} to write or append, default 'w' (optional)
def writeGnuplotFctn(filename, dim, fctn, resolution, mode="w"):
    p = DataVector(dim)
    fout = gzOpen(filename, mode)

    # evaluate 1d function
    if dim == 1:
        for x in xrange(resolution):
                p[0] = float(x) / (resolution - 1)
                pc = fctn(p)
                fout.write("%f %f\n" % (p[0], pc))
    # evaluate 2d function
    elif dim == 2:
        for x in xrange(resolution):
            for y in xrange(resolution):
                p[0] = float(x) / (resolution - 1)
                p[1] = float(y) / (resolution - 1)
                pc = fctn(p)
                fout.write("%f %f %f\n" % (p[0], p[1], pc))
            fout.write("\n")
        fout.write("e\n")
    # can't plot anything else
    else:
        sys.stderr.write("Error! Can't plot grid with dimensionality %d..." % (dim))
    fout.close()
    return

#-------------------------------------------------------------------------------
## @brief Writes coordinates of a grid into a file, suitable for gnuplot.
#
# @param filename Filename to which data is written
# @param grid Grid
def writeGnuplotGrid(filename, grid):
    dim = grid.getDimension()
    gs = grid.getStorage()
    if dim == 2:
        p = DataVector(dim)
        fout = file(filename, "w")
        for i in range(grid.getSize()):
            gp = gs.getCoordinates(gs.getPoint(i), p)
            fout.write("%f %f\n" % (p[0],p[1]))
    # can't plot anything else
    else:
        sys.stderr.write("Error! Can't plot grid with dimensionality %d..." % (dim))
    fout.write("e\n")
    fout.close()
    return

#-------------------------------------------------------------------------------
# read/write alpha
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
## @brief Writes DataVector to arff file. 
# If filename ends with ".gz", file is gzip-compressed.
#
# @param filename Filename of new file
# @param alpha The DataVector
def writeAlphaARFF(filename, alpha):
    fout = gzOpen(filename, "w")
    fout.write("@RELATION \"%s ALPHAFILE\"\n\n" % filename)
    fout.write("@ATTRIBUTE alpha NUMERIC\n")
    
    fout.write("\n@DATA\n")
    
    for i in xrange(len(alpha)):
        fout.write("%1.20f\n" % alpha[i])
    
    fout.close()

#-------------------------------------------------------------------------------
## @brief Reads in onedimensional data from an ARFF file.
#
# @param filename the file's filename that should be red
# @return returns the DataVector
def readAlphaARFF(filename):
    try:
        data = readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        sys.exit(1)

    dv = DataVector(data["data"].getNrows())
    data["data"].getColumn(0, dv)
    return dv

#-------------------------------------------------------------------------------
## @brief Reads in onedimensional data from a delimiter separated data file.
#
# @param filename the file's filename that should be red
# @param delim (optional) separator between columns. Default: whitespaces
# @return returns the DataVector
def readAlphaTrivial(filename, delim = None):
    try:
        data = readDataTrivial(filename, delim, hasclass=False)
    except:
        print ("An error occured while reading " + filename + "!")
        sys.exit(1)

    dv = DataVector(data["data"].getNrows())
    data["data"].getColumn(0, dv)
    return dv

#-------------------------------------------------------------------------------
## @brief Opens and reads the onedimensional data of an ARFF (or plain whitespace-separated data) file.
#
# @param filename filename of the file
# @return the data stored in the file as a set of arrays, DataVector
def readAlpha(filename):
    try:
        if isARFFFile(filename):
            data = readAlphaARFF(filename)
        else:
            data = readAlphaTrivial(filename)
    except Exception, e:
        print ("An error occured while reading " + filename + "!")
        raise e

    return data

#-------------------------------------------------------------------------------
## @brief Serialize a Grid to a file.
# If filename ends with ".gz", file is gzip-compressed.
#
# @param filename Filename of new file
# @param grid The Grid
def writeGrid(filename, grid):
    text = grid.serialize()
    fout = gzOpen(filename, "w")
    fout.write(text)
    fout.close()

#-------------------------------------------------------------------------------
## @brief Unserialize a Grid from a file.
#
# @param filename Filename of file 
# @return Grid
def readGrid(filename):
    fin = gzOpen(filename, "r")
    text = fin.read()
    fin.close()
    
    return Grid.unserialize(text)

#-------------------------------------------------------------------------------
## @brief Write whole checkpoint data to file.
# This writes two files containing grid and coefficient (alpha) vector.
# Optionally, two additional parameters can be specified, influencing the filename.
# The filename has the following form: FILENAME[.aADAPTATION][.fFOLD].{alpha.arff.gz|grid.gz}
#
# @param filename Filename prefix
# @param grid Grid file
# @param alpha Coefficient DataVector
# @param adaption (optional) number of adaptive step for refinement
# @param fold (optional) specifying which fold
def writeCheckpoint(filename, grid, alpha, adaption = None, fold = None):
    adapt_str = ""
    fold_str = ""
    if adaption != None:
        adapt_str = ".a%d" % (adaption)
    if fold != None:
        fold_str = ".f%d" % (fold)
    writeAlphaARFF("%s%s%s.alpha.arff.gz" % (filename, fold_str, adapt_str), alpha)
    writeGrid("%s%s%s.grid.gz" % (filename, fold_str, adapt_str), grid)
    
#-------------------------------------------------------------------------------
# perform split of dataset in num_partitions partitions for n-fold-cv
# split data into folds,
# return ([data1,data2,...,datan], [classes1,classes2,...,classesn])
#-------------------------------------------------------------------------------
def split_n_folds(data, num_partitions, seed=None):
    dim = data["data"].getNcols()
    size = data["data"].getNrows()
    # create permutation
    random.seed(seed)
    seq = range(size)
    random.shuffle(seq)
    # container for new Data and Classes
    dvec = []
    cvec = []
    cv = DataVector(dim)

    size_left = size
    index = 0
    for i in xrange(num_partitions):
        size_fold = size_left/(num_partitions-i)
        dvec.append(DataMatrix(size_fold, dim))
        cvec.append(DataVector(size_fold))
        for rowNum in xrange(size_fold):
            data["data"].getRow(seq[index], cv)
            dvec[i].setRow(rowNum, cv)
            cvec[i][rowNum] = data["classes"][seq[index]]
            index += 1
        size_left = size_left-size_fold
    
    return (dvec, cvec)










# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  Following functions have not yet been updated...
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#-------------------------------------------------------------------------------
## @brief writes statistics
#
# @param filename filename
# @param txt text to write to file
# @param mode writing mode (default "a")
def writeStats(filename, txt, mode = "a"):
    writeLockFile(filename + ".stats.gz", txt, mode)

#-------------------------------------------------------------------------------
## @brief read checkpoint
#
# @param filename filename
def readCheckpoint(filename):
    alpha = readAlphaARFF(filename+".alpha.arff")
    grid = readGrid(filename+".grid")
    
    return grid, alpha


## @brief (Recursively) creates a directory if not yet existant.
#
# @param path Path of directory
# @param verbose Tell what is been done (optional)
def makedir(path, verbose=False):
    if not os.path.isdir(path):
        os.makedirs(path)
        if verbose:
            print "Created directory %s." %(s)
    else:
        if verbose:
            print "Nothing done. Directory %s already existing." %(s)
## @brief write ARFF data
#
# @param data data to write
# @param merge flag to merge data (default False)
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

            if isinstance(dataset["data"], DataMatrix):
                dim = dataset["data"].getNcols()
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
        
        if isinstance(dataset["data"], DataMatrix):
            num_rows = dataset["data"].getNrows()
            for row in xrange(num_rows):
                lout = []
                for column in xrange(dim):
                    lout.append(dataset["data"].get(row,column))
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
## @brief Writes a DataVector object to a file, specified by filename.
#
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
## @brief Writes information that is needed for the normalization of data to a file.
#
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
## @brief Reads information that is needed for the normalization of data from a file.
#
# @param filename a filename
# @return (border, minvals, maxvals, deltavals) @n
# border: offset for normalization @n
# minvals: the minimum value of each attribute @n
# maxvals: the maximum value of each attribute @n
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
## @brief Normalize values of input vectors on the segment [0,1]
#
# @param data Dataset
# @param border Specifies border of the dataset, will be added to the normalized value 
# @param filename Filename of normfile (optional)
# @param minvals Array of normalization boundary min values (one per dimension) (optional)
# @param maxvals Array of normalization boundary max values (one per dimension) (optional)
# @param verbose Provide additional output
def normalize(data, border=0.0, filename=None, minvals=None, maxvals=None, verbose=False):
    # check parameters
    if len(data) == 0:
        raise ValueError, "Wrong or no data."
    if minvals and maxvals:
        if (len(minvals) <> data[0]["data"].getNcols() or 
            len(maxvals) <> data[0]["data"].getNcols()):
            raise ValueError, "Wrong number of min- or max-values."
        lmin = minvals
        lmax = maxvals
    else:
        # init lmin and lmax to first values of first dataset
        lmin = []
        lmax = []
        for dim in range(data[0]["data"].getNcols()):
            lmin.append(data[0]["data"].get(0,dim))
            lmax.append(data[0]["data"].get(0,dim))

        for dataset in data:
            for dim in xrange(dataset["data"].getNcols()):
                cmin = dataset["data"].min(dim)
                lmin[dim] = min(cmin, lmin[dim])

                cmax = dataset["data"].max(dim)
                lmax[dim] = max(cmax, lmax[dim])
    # output
    if verbose:
        print "Dim:", len(lmin)
        print "Boundary:", border
        for d in range(len(lmin)):
            print " [%f,%f]" % (lmin[d], lmax[d])

    # delta values
    ldelta = map(lambda x,y: ((y-x)/(1.0-2.0*border)), lmin, lmax)

    # write normalization data to file:
    if filename:
        if verbose: print "Writing normalization information to", filename
        writeNormfile(filename, border, lmin, lmax)
    
    for dataset in data:
        vec_tmp = DataVector(dataset["data"].getNrows())    
        for dim in xrange(dataset["data"].getNcols()):
            # special handling for the case that all max==min, i.e. all
            # attribute values are equal: set to 0.5
            if ldelta[dim] == 0:
		dataset["data"].getColumn(dim, vec_tmp)
		vec_tmp.setAll(0.5)
                dataset["data"].setColumn(dim, vec_tmp)
            else:
		dataset["data"].getColumn(dim, vec_tmp)
		for j in range(dataset["data"].getNrows()):
		    vec_tmp[j] = (vec_tmp[j]-lmin[dim]) / ldelta[dim] + border
                
                dataset["data"].setColumn(dim, vec_tmp)
    return


#-------------------------------------------------------------------------------
## @brief Divides the class values in two categories
#
# @param data Dataset 
# @param border Classes will be differentiated between greater and less then border 
# @param minborder All classes under the minborder are processed as if they were over border
# @param verbose Provide additional output
def normalizeClasses(data, border=0.0, minborder=-sys.maxint-1, verbose=False):
    if verbose:
        print "Cut-off at", border
    for dataset in data:
        dataset["classes"].partitionClasses(border)
    return

#-------------------------------------------------------------------------------
## @brief Validates Dataset
#
# @param data Dataset 
def checkData(data):
    if len(data) == 0:
    	print("No data loaded. Aborting!")
    	sys.exit(1)

    hasClasses = False
    
    data_len = data[0]["data"].getNcols()
    if data[0].has_key("classes"):
    	hasClasses = True
    
    for dataset in data:
    	if data_len != dataset["data"].getNcols():
    		print("Error! Can't merge data due to different amount of dimensions!")
    		sys.exit(1)
    	if dataset.has_key("classes") != hasClasses:
    		print("Error! Can't merge data, because some files have classification data and some not!")
    		sys.exit(1)
    		
    return



#-------------------------------------------------------------------------------
# perform sequential(!) split of dataset in num_partitions partitions for n-fold-cv
# split data into folds,
# return ([data1,data2,...,datan], [classes1,classes2,...,classesn])
# @param data Dataset 
# @param num_partitions The number of n partitions to split in.
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
def split_DataVector_by_proportion(data, proportion):
    dim = data.getDim()
    size = data.getSize()

    splitpoint = int(min(size-1, round(size*proportion)))
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
## @brief perform stratified split of a data set given by two DataVectors into two DataVectors each
# @param data DataVector with data points to split
# @param classes DataVector with class values to split
# @param proportion split into proportion, (1-proportion)
# @return (data1, data2, classes1, classes2)
def split_DataVectors_by_proportion_stratified(data, 
                                               classes,
                                               proportion):
    proportion = float(proportion)
    dim = data.getDim()
    size = data.getSize()

    splitpoint = int(min(size-1, round(size*proportion)))
    dv1 = DataVector(splitpoint, dim)
    dv2 = DataVector(size-splitpoint, dim)
    cv1 = DataVector(splitpoint, 1)
    cv2 = DataVector(size-splitpoint, 1)
    row = DataVector(1,dim)
    
    # get index ranges for stratification
    index_pos = [i for i in range(size) if classes[i]>=0]
    index_neg = [i for i in range(size) if classes[i]<0]
    index_pos_splitpoint = int(min(size-1, round(len(index_pos)*proportion)))
    index_neg_splitpoint = int(min(size-1, round(len(index_neg)*proportion)))
    
    if (index_pos_splitpoint+index_neg_splitpoint > splitpoint):
        proportion_pos = len(index_pos)*proportion
        proportion_neg = len(index_neg)*proportion
        if proportion_pos-math.floor(proportion_pos) > proportion_neg-math.floor(proportion_neg):
            index_neg_splitpoint -= 1
        else:
            index_pos_splitpoint -= 1

    # copy Data and Classes
    indices1 = index_pos[0:index_pos_splitpoint]+index_neg[0:index_neg_splitpoint]
    indices2 = index_pos[index_pos_splitpoint:]+index_neg[index_neg_splitpoint:]
    i = 0
    for j in indices1:
        data.getRow(j, row)
        dv1.setRow(i, row)
        cv1[i] = classes[j]
        i += 1
    i = 0
    for j in indices2:
        data.getRow(j+splitpoint, row)
        dv2.setRow(i, row)
        cv2[i] = classes[j]
        i += 1

    if len(dv1)+len(dv2) <> len(data): raise Exception("data length doesn't match")
    return (dv1, dv2, cv1, cv2)


def readGridAlpha(fnamegrid, fnamealpha):
    grid = readGrid(fnamegrid)
    grid.getStorage().recalcLeafProperty()
    alpha = readAlpha(fnamealpha)
    return (grid, alpha)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## An array containing all modes and descriptions
CModes = {
    "laplace" : "Classical Laplacian. See OpLaplaceAdaptive",
    "identity" : "Identity matrix, most efficient.",
    "identity_no_level1" : "Identity matrix, most efficient. But do not penalize Level 1",
    "anisotropicpenalty" : "Preferres quadratical supports. See OperationRegularizationDiagonal.",
    "levelsum" : "Sum of the levels, scaled by the gridlevel (usually 2 for adaptive SGs).",
    "isotropicpenalty" : "Energy-norm-like SGs. See OperationRegularizationDiagonal.",
    "rowsum" : "Sum of the rows of classical Laplacian. See OpLaplaceAdaptive",
    "hkmix": "H^k_mix norm; requires parameter Hk",
    "h0hklaplace" : "Pseudo-Laplace with H^k in one, and 'H^0' in the remaining dimensions each; requires paramer Hk"}
    
## base function types
base_types = {
#    "linear" : {"base" : SLinearBase, "b" : SGridOperationB, "t" : test_dataset_linear, "laplace" : SGridOperationLaplace},
#    "modlinear" : {"base" : SLinearModifiedBase, "b" : SGridModOperationB, "t" : test_dataset_modlin},
#    "poly" : {"base" : "SPolyBase", },
              }


## Class Matrix that incorporates settings and actions for applying
# the matrix C and computing the RHS b.
class Matrix:
    def __init__(self, grid, x, l, mode, Hk, base = None):
        self.grid = grid
        self.x = x
        self.l = l
        self.B = createOperationMultipleEval(grid, x)
        #val: base is obviously not used
        #self.base = base
        self.CMode = mode.lower()
        
        if self.CMode == "laplace":
            self.C = createOperationLaplace(grid)
        elif self.CMode == "identity":
            self.C = createOperationIdentity(grid)
        elif self.CMode == "identity_no_level1":
            pass
        elif self.CMode == "anisotropicpenalty":
            self.C = createOperationRegularizationDiagonal(grid, OperationRegularizationDiagonal.ANISOTROPIC_PENALTY, 0)
        elif self.CMode == "levelsum":
            temp = DataVector(grid.getSize())
            # fill temp vector with levelsums
            gridStorage = self.grid.getStorage()
            for i in range(gridStorage.getSize()):
                gp = gridStorage.getPoint(i)
                temp[i] = gp.getLevelSum()
            class Diagop(object):
                def __init__(self, d):
                    self.d = d
                def mult(self, a, res):
                    res.copyFrom(a)
                    res.componentwise_mult(self.d)
            self.C = Diagop(temp)

        elif self.CMode == "isotropicpenalty":
            self.C = createOperationRegularizationDiagonal(grid, OperationRegularizationDiagonal.ISOTROPIC_PENALTY, 0)
        elif self.CMode == "rowsum":
            opL = createOperationLaplace(grid)
            dva = DataVector(grid.getSize())
            dva.setAll(1)
            dres = DataVector(len(dva))
            opL.mult(dva, dres)
            class Diagop(object):
                def __init__(self, d):
                    self.d = d
                def mult(self, a, res):
                    res.copyFrom(a)
                    res.componentwise_mult(self.d)
            self.C = Diagop(dres)
        elif self.CMode == "pseudounit":
            raise Exception("not implemented")
        elif self.CMode == "hkmix":
            self.C = createOperationRegularizationDiagonal(grid, OperationRegularizationDiagonal.HKMIX, Hk)
        elif self.CMode == "h0hklaplace":
            self.C = createOperationRegularizationDiagonal(grid, OperationRegularizationDiagonal.H0HKLAPLACE, Hk)
    
    def generateb(self, y):
        b = DataVector(self.grid.getSize())
        self.B.multTranspose(y, b)
        return b
    
    def ApplyMatrix(self, alpha, result):
        M = self.x.getNrows();
        temp = DataVector(M)
    
        self.B.mult(alpha, temp)
        self.B.multTranspose(temp, result)

        if (self.CMode == "laplace" or
            self.CMode == "hkmix" or 
            self.CMode == "h0hklaplace" or
            self.CMode == "isotropicpenalty" or
            self.CMode == "anisotropicpenalty" or
            self.CMode == "rowsum" or
            self.CMode == "levelsum"):
            temp = DataVector(len(alpha))
            self.C.mult(alpha, temp)
            result.axpy(M*self.l, temp)

        elif self.CMode == "identity":
            result.axpy(M*self.l, alpha)
            
        elif self.CMode == "identity_no_level1":
            result.axpy(M*self.l, alpha)
            # now correct for level 1 again
            gridStorage = self.grid.getStorage()
            gi = GridPoint(gridStorage.getDimension())
            for d in range(gridStorage.getDimension()):
                gi.set(d, 1, 1)
            i = gridStorage.getSequenceNumber(gi)
            result[i] = result[i] - M*self.l*alpha[i]
            
#        elif self.CMode == "levelsum":
#            temp = DataVector(len(alpha))
#            # fill temp vector with levelsums
#            gridStorage = self.grid.getStorage()
#            for i in range(gridStorage.getSize()):
#                gp = gridStorage.getPoint(i)
#                temp[i] = gp.getLevelSum()*alpha[i]
#            result.axpy(M*self.l, temp)
#
#        elif self.CMode == "rowsum":
#            # completely inefficient, but sufficient for test purposes
#            temp = DataVector(len(alpha))
#            ones = DataVector(len(alpha))
#            ones.setAll(1)
#            self.C.mult(ones, temp)
#            for i in range(len(alpha)):
#                temp[i] = temp[i]*alpha[i]
#            result.axpy(M*self.l, temp)
        
        elif self.CMode == "pseudounit":
            raise Exception("not implemented")

        else:
            sys.stderr.write("Error! Mode %s not existant!\n" % (self.CMode))
            sys.exit(1)



# #-------------------------------------------------------------------------------
# # saves/restores grid from file
# # set appropriate modes for automatic destructions
# #-------------------------------------------------------------------------------
# 
# ## @brief
# def restoreGrid(text):
#     return Grid.unserialize(text)
# 
# 
# def saveGrid(grid):
#     return grid.serialize()
# 
##-------------------------------------------------------------------------------
### @brief Converts a python "dataset"-structure into two objects of type DataVector (X,Y).
##
#def createDataVectorFromDataset(dataset):
#    dim = len(dataset["data"])
#    entries = len(dataset["data"][0])
#    data = dataset["data"]
#    
#    x = DataVector(entries, dim)
#    y = None
#
#    for d in xrange(dim):    
#        for i in xrange(entries):
#            x[i*dim + d] = data[d][i]
#            
#    if dataset.has_key("classes"):
#        classes = dataset["classes"]
#        y = DataVector(entries)
#        for i in xrange(entries):
#            y[i] = classes[i]
#            
#    return (x,y)
#
##-------------------------------------------------------------------------------
### @brief Converts one or two (data + optionally classes) objects of type DataVector
## to a python "dataset"-structure.
#def createDatasetFromDataVector(data, classes=None):
#    dataset = {}
#    dataset['data'] = []
#    for j in range(data.getDim()):
#        column = [data.get(i,j) for i in range(data.getSize())]
#        dataset['data'].append(column)
#
#    if classes:
#        dataset['classes'] = [classes[i] for i in range(len(classes))]
#    
#    return dataset

#-------------------------------------------------------------------------------
## @brief Read in an arbitrary grid
