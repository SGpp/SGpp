#!/usr/bin/env python
# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# This python script installs the pysgpp package using setuptools. The
# script creates a pysgpp lib in the site-packages of
# python. Furthermore, it collects all the relevant python code,
# located in each module under the path <module name>/python and
# copies it to the site-package folder of pysgpp using the following
# scheme: <path to site-packages>/pysgpp-<unique
# key>/extensions/<module name>/<copy of python code of correponding
# module>.

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
with open(initFile, "w") as f:
    f.write("""
# add current directory to PYTHONPATH such that pysgpp_swig can be imported
import os
import sys
sys.path.append(os.path.dirname(__file__))

# import pysgpp_swig and extensions
from pysgpp_swig import *
import pysgpp.extensions
""")

if len(moduleFolders) > 0:
    # create __init__.py file which imports all the extensions
    initFile = os.path.join("__init__.py")
    with open(initFile, "w") as f:
        for moduleFolder, _ in pythonModuleFolders:
            f.write("import %s\n" % moduleFolder)

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
