# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import os
import sys

# get all subdirs of path, required by CheckJNI
def getSubdirs(path):
    pathlist = []
    for f in os.listdir(path):
        if os.path.isdir(os.path.join(path, f)):
            pathlist.append(os.path.join(path, f))
    return pathlist

# Custom test for executables used during configuration
def CheckExec(context, cmd):
    context.Message('Checking for %s...' % (cmd))
    ret = context.env.WhereIs(cmd)
    if ret == None:
        ret = ''
    context.Result(ret)
    return ret

def CheckFlag(context, flagString):
    
    checkSrcFile = """
int main(int argc, char **argv) {
    /**
     * Does nothing, just a test whether a compilation with a specific flag is possible.
     */
}
    """
    context.Message('Checking for flag "' + flagString + '" ... ')
    oldFlags = context.env['CPPFLAGS']
    context.env.AppendUnique(CPPFLAGS = flagString)
    result = context.TryCompile(checkSrcFile, ".cpp")
    context.Result(result)
    context.env.Replace(CPPFLAGS = oldFlags)
    return result

# Check for jni header file
# if found, additionally add all subdirs to CPPPATH (platform dependent files)
def CheckJNI(context):
    print "Trying to locate jni.h..."
    # message if JNI_CPPINCLUDE not set
    if not os.environ.get('JNI_CPPINCLUDE'):
        print "... JNI_CPPINCLUDE not set"
    # check for JAVA_HOME first
    if os.environ.get('JAVA_HOME'):
        pname = os.path.join(os.environ.get('JAVA_HOME'), 'include')
        if os.path.exists(os.path.join(pname, 'jni.h')):
            context.env.Append(CPPPATH=[pname] + getSubdirs(pname))
            res = "... found in " + pname
            context.Result(res)
            return res
        else:
            print "... not found in $JAVA_HOME/include"
    else:
        print "... JAVA_HOME not set"

    # not found, try guessing:
    # look, where java and javac are located:
    # include/ directory might be 1 or 2 dirs below
    print "... trying to guess"
    for f in ['java', 'javacc']:
        fdir = context.env.WhereIs(f)
        if not fdir:
            continue
        # os.path.realpath to resolve links
        basedir = os.path.dirname(os.path.realpath(fdir))
        for subdir in ['..', os.path.join('..', '..')]:
            pname = os.path.join(basedir, subdir, 'include')
            if os.path.exists(os.path.join(pname, 'jni.h')):
                context.env.Append(CPPPATH=[pname] + getSubdirs(pname))
                res = "... found in " + pname
                context.Result(res)
                return res
            
    context.Result('... nothing found!')
    return 0
