import glob
import os
import sys
import subprocess
import shlex
from matplotlib.compat.subprocess import CalledProcessError

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