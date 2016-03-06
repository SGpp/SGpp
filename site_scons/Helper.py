# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import glob
import os
import re
import sys

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
    languageSupport = []

    for modulePath in modulePaths:
        module = modulePath[:-len(suffix)]
        module = module[len(path):]
        if module in ignoreFolders:
            continue
        if module in ['jsgpp', 'pysgpp']:
          languageSupport.append(module)
          continue
        modules.append(module)
    return modules, languageSupport

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
        with open(env['CMD_LOGFILE'], 'a') as logFile:
            logFile.write('%s\n' % s)

#creates a Doxyfile containing proper module paths based on Doxyfile_template
def prepareDoxyfile(modules):
    '''Create Doxyfile(s) and overview-pages
    @param modules list of modules'''

    # create Doxyfile
    with open('Doxyfile_template', 'r') as doxyFileTemplate:
        with open('Doxyfile', 'w') as doxyFile:
            inputPaths = 'INPUT ='
            excludePaths = 'EXCLUDE ='
            examplePaths = 'EXAMPLE_PATH ='
            imagePaths = 'IMAGE_PATH ='

            for moduleName in modules:
                inputPath = moduleName + '/'
                examplePath = moduleName + '/examples'
                testPath = moduleName + '/tests'
                imagePath = moduleName + '/doc/doxygen/images'

                #print os.path.join(os.getcwd(),inputPath)
                if os.path.exists(os.path.join(os.getcwd(), inputPath)):
                    inputPaths += " " + inputPath
                if os.path.exists(os.path.join(os.getcwd(), examplePath)):
                    examplePaths += " " + examplePath
                    excludePaths += " " + examplePath
                if os.path.exists(os.path.join(os.getcwd(), testPath)):
                    excludePaths += " " + testPath
                if os.path.exists(os.path.join(os.getcwd(), imagePath)):
                    imagePaths += " " + imagePath

            for line in doxyFileTemplate.readlines():
                if re.match(r'INPUT  .*', line):
                    doxyFile.write(inputPaths + "\n")
                elif re.match(r'EXCLUDE  .*', line):
                    doxyFile.write(excludePaths + "\n")
                elif re.match(r'EXAMPLE_PATH  .*', line):
                    doxyFile.write(examplePaths + "\n")
                elif re.match(r'IMAGE_PATH  .*', line):
                    doxyFile.write(imagePaths + "\n")
                else:
                    doxyFile.write(line)

    # create example menu page
    with open('base/doc/doxygen/examples.doxy', 'w') as examplesFile:
        examplesFile.write('''/**
@page examples Examples

This is a collection of examples from all modules.

If you're new to SG++ or want to try out quickly,
read the @ref code_examples_tutorial first.

To add new examples to the documentation,
go to the respective folder MODULE_NAME/doc/doxygen/ and
add a new example file code_examples_NAME.doxy with doxygen-internal
name code_examples_NAME.

Note that SCons automatically compiles (but not runs)
all C++ examples on each run.
For this to work, the examples must lie in the directories of the form
\c /path/to/SGpp/trunk/MODULE_NAME/examples.

''')

        modules.sort()
        tutorial = 'code_examples_tutorial'
        
        for moduleName in modules:
            examplesFile.write('<h2>Module {}</h2>\n'.format(moduleName))
            subpages = glob.glob(os.path.join(
                moduleName, 'doc', 'doxygen', 'code_examples_*.doxy'))
            subpages = [os.path.split(path)[-1][:-5]
                        for path in glob.glob(os.path.join(
                            moduleName, 'doc', 'doxygen',
                            'code_examples_*.doxy'))]
            subpages.sort()
            if tutorial in subpages:
                del subpages[subpages.index(tutorial)]
                subpages = [tutorial] + subpages
            
            for subpage in subpages:
                examplesFile.write('- @subpage {}\n'.format(subpage))

        examplesFile.write('**/\n')

    # create module page
    with open('base/doc/doxygen/modules.doxy', 'w') as modulesFile:
        with open('base/doc/doxygen/modules.stub0', 'r') as stubFile:
            modulesFile.write(stubFile.read())

        for moduleName in modules:
            for subpage in glob.glob(os.path.join(moduleName, 'doc', 'doxygen', 'module_*.doxy')):
                modulesFile.write('- @subpage ' + os.path.splitext(os.path.split(subpage)[-1])[0] + '\n')

        with open('base/doc/doxygen/modules.stub1', 'r') as stubFile:
            modulesFile.write(stubFile.read())


def flatDependencyGraph(dependencies, acc):
    for dependency in dependencies[::-1]:
        if dependency not in acc:
            acc = [dependency] + acc
    return acc



# On win32, the command lines are limited to a ridiculously short length
# (1000 chars). However, compiler/linker command lines easily exceed that
# length. The following is a fix for that.
# It has to be enabled with "env['SPAWN'] = win32_spawn".
# (see https://bitbucket.org/scons/scons/wiki/LongCmdLinesOnWin32)
def set_win32_spawn(env):
    import win32file
    import win32event
    import win32process
    import win32security

    def win32_spawn(sh, escape, cmd, args, spawnenv):
        for var in spawnenv:
            spawnenv[var] = spawnenv[var].encode('ascii', 'replace')

        sAttrs = win32security.SECURITY_ATTRIBUTES()
        StartupInfo = win32process.STARTUPINFO()
        newargs = ' '.join(map(escape, args[1:]))
        cmdline = cmd + " " + newargs

        # check for any special operating system commands
        if cmd == 'del':
            for arg in args[1:]:
                win32file.DeleteFile(arg)
            exit_code = 0
        else:
            # otherwise execute the command.
            try:
              hProcess, hThread, dwPid, dwTid = win32process.CreateProcess(
                None, cmdline, None, None, 1, 0, spawnenv, None, StartupInfo)
            except:
              import win32api
              error_code = win32api.GetLastError()
              raise RuntimeError("Could not execute the following " + 
                  "command line (error code {}): {}".format(
                    error_code, cmdline))
            win32event.WaitForSingleObject(hProcess, win32event.INFINITE)
            exit_code = win32process.GetExitCodeProcess(hProcess)
            win32file.CloseHandle(hProcess)
            win32file.CloseHandle(hThread)
        return exit_code
    
    env['SPAWN'] = win32_spawn
