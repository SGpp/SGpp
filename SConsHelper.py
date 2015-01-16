import glob
import os
import sys
import subprocess
import shlex
from matplotlib.compat.subprocess import CalledProcessError

def setupSelectedModules(env, allModules, allModulesNames, allLanguageWrappersList):
    modulesToBuild = []
    modulesToBuildFolder = []
    # if neither module nor support language set, do all
    anySet = False
    for entry in allModulesNames + allLanguageWrappersList:
        if env[entry]:
            anySet = True
      
    # for clean enable everything:
    if env.GetOption('clean') or env['SG_ALL'] or not anySet:
        env['SG_ALL'] = True
        for entry in allModulesNames + allLanguageWrappersList:
            modulesToBuild += [entry]
            env[entry] = True
        for module in allModules:
            if env['SG_' + module.upper()]:
                modulesToBuildFolder += [module]
    return modulesToBuild, modulesToBuildFolder

def setupCompilerDefines(env, modulesToBuild):
    # add C++ defines for all modules
    cppdefines = []
    for module in modulesToBuild:
        if env[module]:
            cppdefines.append(module)
    print cppdefines
    env.Append(CPPDEFINES=cppdefines)

# get all files in a folder matching "SConscript*"
# path has to end with "/"
def getModules(path):
    if path[-1] != '/':
        path += '/'
    
    suffix = '/SConscript'
    searchString = path + '*' + suffix
    modulePaths = glob.glob(searchString)
    modules = []
    for modulePath in modulePaths:
        module = modulePath[:-len(suffix)]
        module = module[len(path):]        
        modules.append(module)
    return modules

# Definition of flags / command line parameters for SCons
#########################################################################

def multiParamConverter(s):
    print s
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