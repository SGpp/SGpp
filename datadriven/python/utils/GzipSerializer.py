# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import gzip, re

##The class provides the method, used by subclasses to work with gzip-compressed
# files 
class GzipSerializer(object):
    
    ## Opens a file. If the file ends with ".gz", automatically gzip compression
    # is used for the file. Returns the filedescriptor
    # @param filename
    # @param mode access mode, default: "r" for read only
    # @return file descriptor
    def gzOpen(self, filename, mode="r"):
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
    
    def serialize(self, targetObj, stream):
        text = str(targetObj)
        stream.write(text)
        
    def serializeToFile(self, targetObj, filename):
        stream = self.gzOpen(filename, "w")
        try:
            self.serialize(targetObj, stream)
        finally:
            stream.close()
            
    def deserialize(self, stream):
        text = stream.read()
        return text
    
    def deserializeFromFile(self, filename):
        result = None
        stream = self.gzOpen(filename, "r")
        try:
            result = self.deserialize(stream)
        finally:
            stream.close()
        return result