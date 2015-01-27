# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import glob
import os
import sys
import subprocess
import shlex
import fileinput
import re

# get all folders containing an "SConscript*" file
# path has to end with "/"
def getModules(ignoreFolders):
    
#     if path[-1] != '/':
#         path += '/'
    #path = '#/'
    path = ''
    suffix = '/SConscript'
    searchString = path + '*' + suffix
    modulePaths = glob.glob(searchString)
    modules = []
    for modulePath in modulePaths:
        module = modulePath[:-len(suffix)]
        module = module[len(path):]
        if module in ignoreFolders:
            continue        
        modules.append(module)
    return modules

# Definition of flags / command line parameters for SCons
#########################################################################

def multiParamConverter(s):
    return s.split(',')

# detour compiler output
def print_cmd_line(s, target, src, env):
    if env['VERBOSE']:
        sys.stdout.write(u'%s\n' % s)
    else:
        sys.stdout.write(u'.')
        sys.stdout.flush()
    if env['CMD_LOGFILE']:
        open(env['CMD_LOGFILE'], 'a').write('%s\n' % s);

#creates a Doxyfile containing proper module paths based on Doxyfile_template
def prepareDoxyfile(modules):
    
    doxyFileTemplate = open('Doxyfile_template','r')
    doxyFile = open('Doxyfile', 'w')

    inputPaths = 'INPUT = bin/ doc/ '
    examplePaths = 'EXAMPLE_PATH = '
    imagePaths = 'IMAGE_PATH = '

    for moduleName in modules:
        inputPaths += moduleName + '/ '
        examplePaths += moduleName + '/examples '
        imagePaths += moduleName + '/doc/doxygen/images '
	
    for line in doxyFileTemplate.readlines():
        if re.search(r'INPUT  .*', line):	   
	    doxyFile.write(re.sub(r'INPUT.*', inputPaths, line))
        elif re.search(r'EXAMPLE_PATH  .*', line):
	    doxyFile.write(re.sub(r'EXAMPLE_PATH.*', examplePaths, line))
        elif re.search(r'IMAGE_PATH  .*', line):
	    doxyFile.write(re.sub(r'IMAGE_PATH.*', imagePaths, line))
        else:
            doxyFile.write(line)

    doxyFile.close()