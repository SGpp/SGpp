﻿#!/usr/bin/env python

import os
import shutil
from setuptools import setup

# path to pysgpp lib
libpath = os.path.join("lib", "pysgpp")

# list of all available modules -> all folders in the root directory
moduleFolders = [filename for filename in os.listdir(".")
                 if os.path.isdir(filename)]

# select the ones which have a python extension
pythonModuleFolders = [(moduleFolder, os.path.join(moduleFolder, "python"))
                       for moduleFolder in moduleFolders
                       if os.path.exists(os.path.join(moduleFolder, "python"))]

# create the data file list such that it can be used by setuptools
dataFiles = []
for moduleFolder, srcdir in pythonModuleFolders:
    basepath = os.path.join("pysgpp", "extensions", moduleFolder)
    for root, dirs, files in os.walk(srcdir):
        if '.svn' in dirs:
            dirs.remove('.svn')

        dataFiles += [(root.replace(srcdir, basepath),
                       [os.path.join(root, f) for f in files])]

# write init file for pysgpp
initFile = os.path.join(libpath, "__init__.py")
fd = open(initFile, "w")
fd.write("from pysgpp_swig import *%s" % os.linesep)
if len(moduleFolders) > 0:
    fd.write("import extensions")
fd.close()

if len(moduleFolders) > 0:
    # create __init__.py file which imports all the extensions
    initFile = os.path.join("__init__.py")
    fd = open(initFile, "w")
    for moduleFolder, _ in pythonModuleFolders:
        fd.write("import %s%s" % (moduleFolder, os.linesep))
    fd.close()

    dataFiles += [(os.path.join("pysgpp", "extensions"), [initFile])]

# if the current system is windows we need to rename the dll to pyd
dllLibs = [filename for filename in os.listdir(libpath)
          if filename.endswith("dll")]

for dllLib in dllLibs:
    pydLib = "%s.pyd" % os.path.splitext(dllLib)[0]
    shutil.copy2(os.path.join('lib', 'pysgpp', dllLib),
                 os.path.join('lib', 'pysgpp', pydLib))

# setup pysgpp
setup(name='pysgpp',
      version="1.0.0",
      url='sgpp.sparsegrids.org',
      author="Fabian Franzelin",
      description='',
      license='',
      long_description="README",
      platforms='any',
      zip_safe=False,
      package_dir={'': 'lib'},
      packages=['pysgpp'],
      package_data={'pysgpp': ['*.so', '*.lib', '*.pyd']},
      data_files=dataFiles
      )

# cleanup
if len(moduleFolders) > 0 and os.path.exists(initFile):
    os.remove(initFile)
