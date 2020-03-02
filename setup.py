#!/usr/bin/env python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# This python script installs the pysgpp package using setuptools. The
# script creates a pysgpp lib in the site-packages of
# python. Furthermore, it collects all the relevant python code,
# located in each module under the path <module name>/python and
# copies it to the site-package folder of pysgpp under the extensions namespace
# pysgpp.extensions

import os
import shutil
from setuptools import setup, find_packages


try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False
except ImportError:
    bdist_wheel = None

# path to pysgpp lib
libpath = os.path.join("lib", "pysgpp")
extensionspath = os.path.join(libpath, "extensions")

# list of all available modules -> all folders in the root directory
moduleFolders = [filename for filename in os.listdir(".")
                 if os.path.isdir(filename)]

# select the ones which have a python extension
pythonModuleFolders = [(moduleFolder, os.path.join(moduleFolder, "python"))
                       for moduleFolder in moduleFolders
                       if os.path.exists(os.path.join(moduleFolder, "python"))]

try:
    os.mkdir(extensionspath)
except FileExistsError as e:
    pass

# create list of extension scripts
extFiles = []
for moduleFolder, srcdir in pythonModuleFolders:
    basepath = os.path.join(extensionspath, moduleFolder)
    try:
        os.mkdir(basepath)
    except FileExistsError as e:
        pass
    for root, dirs, files in os.walk(srcdir):
        if '.git' in dirs:
            dirs.remove('.git')

        extFiles += [os.path.join(root, f) for f in files if ".py" in f]

##
# copy extension python files to new layout
# pysgpp
# --extensions
# ----modulename
# ------*.py
##

for f in extFiles:
    dest = os.path.join(extensionspath,f)
    dest = dest.replace(os.sep + "python", "")
    try:
        shutil.copy2(f, dest)
    except FileNotFoundError as e:
        os.mkdir(os.path.dirname(dest))
        shutil.copy2(f, dest)
    except shutil.SameFileError as e:
        pass
    

# write init file for pysgpp
initFile = os.path.join(libpath, "__init__.py")
with open(initFile, "w") as f:
    f.write("""
# add current directory to PYTHONPATH such that pysgpp_swig can be imported
import os
import sys
sys.path.append(os.path.dirname(__file__))

# import pysgpp_swig and extensions
from .pysgpp_swig import *
from . import extensions
""")

if len(moduleFolders) > 0:
    # create __init__.py file which imports all the extensions
    initFile = os.path.join("__init__.py")
    with open(initFile, "w") as f:
        for moduleFolder, _ in pythonModuleFolders:
            f.write("from . import %s\n" % moduleFolder)

    try:
        shutil.copy2(initFile, os.path.join(extensionspath, initFile))
    except shutil.SameFileError as e:
        pass
    

# if the current system is windows we need to rename the dll to pyd
dllLibs = [filename for filename in os.listdir(libpath)
          if filename.endswith("dll")]

for dllLib in dllLibs:
    pydLib = "%s.pyd" % os.path.splitext(dllLib)[0]
    shutil.copy2(os.path.join('lib', 'pysgpp', dllLib),
                 os.path.join('lib', 'pysgpp', pydLib))

# setup pysgpp
setup(name='pysgpp',
      version="0.0.0",
      url='https://github.com/SGpp/SGpp',
      author="Dirk.Pflueger@ipvs.uni-stuttgart.de",
      description='''The sparse grids toolkit SG++
 SG++ is a collection of numerical algorithms for sparse grids. It
 contains modules for interpolation, quadrature, data mining
 (regression, classification, clustering), optimization, PDEs, and
 more. SG++ implements algorithms for spatially adaptive grids and
 also provides a module for the combination technique. Many of the
 implemented algorithms are also available as a high-performance
 version, often orders of magnitude faster than standard
 implementations.''',
      license='BSD-style license',
      long_description="README",
      zip_safe=False,
      package_dir={'': 'lib'},
      packages=find_packages(where='lib', include=['pysgpp', 'pysgpp.extensions*']),
      package_data={'pysgpp': ['_pysgpp_swig.so', '*.lib', '*.pyd']},
      )

# cleanup
if len(moduleFolders) > 0 and os.path.exists(initFile):
    os.remove(initFile)
    shutil.rmtree(extensionspath, ignore_errors=True)
