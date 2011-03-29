# file: addnamespaceusing.py
from os.path import *
import re

def main():
    dirdump = '__dirdump.txt'
    namespace = "datadriven"
    namespace_filename = namespace + 'files.txt'
    b = open(namespace_filename)
    namespace_files = b.readlines()
    b.close()
    classnames = {}
    for line in namespace_files:
        if isfile('src/sgpp/' + line.strip()):
            filename = basename(line.strip())
            classname = filename[:filename.index('.')]
            classnames[classname] = True
    classnames = set(classnames.keys())
    #print classnames

    f = open(dirdump)
    for line in f.readlines():
        filename = 'src/sgpp/' + line.strip()
        if isfile(filename) and line not in namespace_files:
            codefile = open(filename)
            code = codefile.read()
            codefile.close()
            codefile = open(filename)
            codelines = codefile.readlines()
            codewords = re.split('\W+', code)
            codefile.close()
            # if file contains classes from namespace package
            if classnames.intersection(codewords):
                s = ""
                using_added = False
                
                for i in xrange(len(codelines)-1, -1, -1):
                    if "using namespace sg::%s;" % namespace in codelines[i]:
                        break
                    if "#include " in codelines[i] and not using_added:
                        s = codelines[i] + ( "using namespace sg::%s;\n" % namespace ) + s
                        using_added = True
                    else:
                        #if "sg::" in codelines[i]:
                        #    print "%s:%d:%s" % (filename, i, codelines[i].strip())
                        s = codelines[i].replace("sg::",'') + s
                        
                if using_added:
                    print filename, "changed"
                    codefile = open(filename,'w')
                    codefile.write(s)
                    codefile.close()
                        
            
        
        
if __name__ == '__main__':
    main()
            
