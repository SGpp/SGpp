


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
import sys, os
sys.path.append(os.path.abspath(os.path.join(sys.path[0], "..", "..")))
from bin.tools import *

if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("-i", "--infile", action="append", type="string", dest="infiles", 
              help="Specifies the inputfile to normalize. If multiple inputfiles are given they are normalized together, not separately.")
    parser.add_option("-t", "--types", action="append", type="string", dest="types", 
              help="Specifies the type of the inputfiles. Availible types are \"arrf\" and \"simple\". Can be specified for each inputfile.")
    parser.add_option("-m", "--merge", action="store_true", default=False, dest="merge", 
              help="If this is enabled, then all inputfiles are normalized and merged to one outputfile.")
    parser.add_option("-n", "--nonormalization", action="store_true", default=False, dest="nonormalization", 
              help="If this option is set, neither data (-b) nor classes (-c) are normalized, but all data is just being converted.")
    parser.add_option("--nodatanormalization", action="store_true", default=False, dest="nodatanormalization", 
              help="If this option is set, data (-b) is not normalized, but just being converted.")
    parser.add_option("--noclassnormalization", action="store_true", default=False, dest="noclassnormalization", 
              help="If this option is set, classes (-c) are not normalized, but just being converted.")
    parser.add_option("-b", "--border", action="store", type="float", dest="border", 
              help="Specifies the border for the Dataset. If not set 0.0 is used.", default=0.0)
    parser.add_option("-c", "--class", action="store", type="float", dest="c_border",metavar="BORDER", 
              help="Specifies the classification border on which the classification data is put into different classes. If not set classvalue>=0.5 is used.", default=0.5)
    parser.add_option("--class_min", action="store", type="float", dest="c_border_min", default=-sys.maxint-1, metavar="BORDERMIN", 
              help="Specifies a second classification border on which the classification data is put into different classes (additionally check for classvalue<=class_min). If not set, it is ignored.")
    parser.add_option("-C", "--noclasses", action="store_true", default=False, dest="noclasses", 
              help="If this is enabled, then inputfiles have no classes.")
    parser.add_option("-o", "--output", action="append", type="string", dest="outfiles", 
              help="Specifies the output file. Can be specified for each inputfile. If not applicated, .arff.gz is appended to each file.")
    parser.add_option("--delimiter", action="store", type="string", default="", dest="delimiter", 
              help="The delimiter separating the columns for the simple-format. Default: Split for whitespaces.")
    parser.add_option("--normfile", action="store", type="string", dest="normfile", default=None, 
              help="Write normalization information to file, so that further data could be normalized")
    parser.add_option("--maple", action="store_true", default=False, dest="maple", 
              help="If enabled, write Maple-readable format.")
    parser.add_option("--hasnoclass", action="store_false", default=True, dest="hasnoclass", 
              help="Set, if data does not contain class attribute.")
    parser.add_option("--min", action="append", type="float", default=None, dest="min", 
              help="Set min value of feature space (normalization boundary). One entry per dimension and in the correct order.")
    parser.add_option("--max", action="append", type="float", default=None, dest="max", 
              help="Set max value of feature space (normalization boundary). One entry per dimension and in the correct order.")
    parser.add_option("--omitt", action="append", type="int", default=None, dest="omitt", 
              help="Can be set multiple times to omitt certain attributes (0,1,2,...)")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", 
              help="Provides extra output")
    (options,args)=parser.parse_args()
    
    if options.infiles == None:
        print("No inputfiles. Aborting...")
        print
        parser.parse_args(["-h"])
    
    if len(args) != 0:
        print("Warning: There were command-line args. Maybe forgotten to put an -i in front of filename?")
    
    if options.outfiles != None and len(options.outfiles) > len(options.infiles):
        print("Warning: There are more outputfiles than inputfiles. Excess arguments are ignored.")
        options.outfiles = options.outfiles[:len(options.infiles)]
    
    if options.types != None and len(options.types) > len(options.infiles):
        print("Warning: There are more types specified than inputfiles. Excess arguments are ignored.")
        options.types = options.types[:len(options.infiles)]
    
    if options.outfiles == None:
        options.outfiles = []
    
    if options.types == None:
        options.types = []
    
    # Now all arrays have been created and are only as lang as infile[]
    
    for i in xrange(len(options.infiles)):
        if i >= len(options.outfiles):
            # append ".arff.gz", but strip trailing ".gz"
            if options.infiles[i].endswith(".gz"):
                options.outfiles.append(options.infiles[i][:-3] + ".arff.gz")
            else:
                options.outfiles.append(options.infiles[i] + ".arff.gz")
    
        if i >= len(options.types):
            filename = options.infiles[i].lower()
            if filename.endswith("arff") or filename.endswith("arff.gz"):
                options.types.append("arff")
            else:
                options.types.append("simple")
        else:
            options.types[i] = options.types[i].lower()
            if not (options.types[i] == "arff" or options.types[i] == "simple"):
                print("Filetype " + options.types[i] + " is unknown. Aborting ...")
                sys.exit(1)
    
    # Gathered all information:
    # in options.infiles[] there are the inputfiles,
    # in options.outfiles[] the output files
    # in options.types[] all types
    
    data = []
    
    for i in xrange(len(options.infiles)):
        try:
#            if options.types[i] == "arff":
#                data.append(readDataARFF(options.infiles[i]))
#                data[i]["filename"] = options.outfiles[i]
#            elif options.types[i] == "simple":
#                data.append(readDataTrivial(options.infiles[i], delim=options.delimiter, hasclass=options.hasnoclass))
#                data[i]["filename"] = options.outfiles[i]
            data.append(readData(options.infiles[i]))
            data[i]["filename"] = options.outfiles[i]
        except Exception, e:
            print("Error while reading "  + options.infiles[i] +"! Aborting...");
            print e
            sys.exit(1)
    
    # sanity check                
    checkData(data)
    
    # omitt certain attributes?
    if options.omitt:
        print "Omitting out of function!"
        dim = data[0]["data"].getNcols()
        attribs = range(dim)
        for a in options.omitt:
            attribs.remove(a)
        if options.verbose: print "Omitting attribute(s)", options.omitt, " taking", attribs
        for dataset in data:
            dataset["data"] = [dataset["data"][i] for i in attribs]
    
    if not options.nonormalization:
        if not options.nodatanormalization:
            if options.verbose: print "Normalizing data..."
            normalize(data, options.border, options.normfile, options.min, options.max, options.verbose)
        if not options.noclassnormalization:
            if options.verbose: print "Normalizing classes..."
            normalizeClasses(data, options.c_border, options.c_border_min, options.verbose)
    if not options.maple:
        writeDataARFF(data, options.merge)
    else:
        writeDataMaple(data, options.merge)