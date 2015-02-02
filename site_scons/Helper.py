# Copyright (C) 2008-today The SG++ Project
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
    '''Create Doxyfile(s) and overview-pages
    @param modules list of modules'''
    
    # create Doxyfile
    doxyFileTemplate = open('Doxyfile_template','r')
    doxyFile = open('Doxyfile', 'w')

    inputPaths = 'INPUT = bin/ doc/'
    examplePaths = 'EXAMPLE_PATH ='
    imagePaths = 'IMAGE_PATH ='

    for moduleName in modules:
        
        inputPath = moduleName + '/'
        examplePath = moduleName + '/examples'
        imagePath = moduleName + '/doc/doxygen/images'
        
        print os.path.join(os.getcwd(),inputPath)
        if os.path.exists(os.path.join(os.getcwd(),inputPath)):
            inputPaths += " " + inputPath
        if os.path.exists(os.path.join(os.getcwd(),examplePath)):
            examplePaths += " " + examplePath
        if os.path.exists(os.path.join(os.getcwd(),imagePath)):
            imagePaths += " " + imagePath
	
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
    
    # create example menu page
    examplesFile = open('base/doc/doxygen/examples.doxy', 'w')
    examplesFile.write('/**\n')
    examplesFile.write('@page examples Examples\n\n')
    examplesFile.write('This is a collection of examples from all modules.\n')
    examplesFile.write('To add new examples, go to the respective folder module/doc/doxygen/\n')
    examplesFile.write('and add a new example file code_examples_NAME.doxy with doxygen-internal\n')
    examplesFile.write('name code_examples_NAME.\n\n')
    for moduleName in modules:
        for subpage in glob.glob(os.path.join(moduleName, 'doc', 'doxygen', 'code_examples_*.doxy')):
            examplesFile.write('- @subpage ' + (os.path.split(subpage)[-1])[:-5] + '\n')
    examplesFile.write('**/\n')

    # create module page
    modulesFile = open('base/doc/doxygen/modules.doxy', 'w')
    modulesFile.write(open('base/doc/doxygen/modules.stub0', 'r').read())
    for moduleName in modules:
        for subpage in glob.glob(os.path.join(moduleName, 'doc', 'doxygen', 'module_*.doxy')):
            modulesFile.write('- @subpage ' + (os.path.split(subpage)[-1])[:-5] + '\n')
    modulesFile.write(open('base/doc/doxygen/modules.stub1', 'r').read())
    modulesFile.close()
