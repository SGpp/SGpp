#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################


import os
import distutils.sysconfig

opts = Options('custom.py')

opts.Add('CPPFLAGS','Set additional Flags','')
opts.Add('LINKFLAGS','Set additional Linker-flags','')

opts.Add('MARCH','Set processor specific MARCH', None)

opts.Add('ICC', 'Uses Intels Optimizing Compiler', False)
opts.Add('OMP', 'Use OpenMP parallelisation', False)
opts.Add('INTELHOME', 'Intel Compiler Home Dir', '')

opts.Add('JSGPP', 'Build jsgpp if set to True', False)
opts.Add('JNI_CPPPATH', 'Path to JNI includes', None)
opts.Add('JNI_OS', 'JNI os path', None)

env = Environment(options = opts, ENV = os.environ)

env.Append(CPPFLAGS=['-pthread'])
env.Append(CPPFLAGS=['-Wall', '-ansi', '-Wno-long-long']) #, '-pedantic'
env.Append(LINKFLAGS=['-pthread'])

####### enable omp support #######
if not env['ICC']:
    if env['OMP']:
        env.Append(CPPFLAGS=['-fopenmp'])
        env.Append(CPPDEFINES=['USEOMP'])
        env.Append(LINKFLAGS=['-fopenmp'])

####### Sets icc tool chain #######
if env['ICC']:
#    env.Tool(env['INTELHOME'] + 'icpc')
    env['CC'] = (env['INTELHOME'] + 'icc')
    env['LINK'] = (env['INTELHOME'] + 'icpc')
    env['CXX'] = (env['INTELHOME'] + 'icpc')
    
if env['ICC']:
    if env['OMP']:
        env.Append(CPPFLAGS=['-openmp'])
        env.Append(CPPDEFINES=['USEOMP'])
        env.Append(LINKFLAGS=['-openmp'])    

if not env['ICC'] and env.has_key('MARCH'):
	env.Append(CPPFLAGS=('-march=' + env['MARCH']))


env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])


if not env.GetOption('clean'):	
    config = env.Configure()
	
    if env['ICC']:
        if not config.CheckLib('guide'):
            Exit(1)

    if env['ICC']:
        if not config.CheckLib('svml'):
            print "SVML should be available when using intelc. Consider runnning scons --config=force!"

    if not config.CheckLibWithHeader('m', 'math.h', 'c++'):
        Exit(1)
        
    if not config.CheckCHeader('Python.h'):
    	print "Python not found. Check path to Python include files."
    	Exit(1)

    env = config.Finish()

Export('env')

SConscript('src/sgpp/SConscript', build_dir='tmp/build_sg', duplicate=0)
SConscript('src/pysgpp/SConscript', build_dir='tmp/build_pysgpp', duplicate=0)

if env['JSGPP']:
    SConscript('src/jsgpp/SConscript', build_dir='tmp/build_jsgpp', duplicate=0)

SConscript('tests/SConscript')

cpy = []
cpy += Command("#lib/pysgpp/_pysgpp.so", "#/tmp/build_pysgpp/_pysgpp.so", Copy("$TARGET", "$SOURCE"))
cpy += Command("#lib/pysgpp/pysgpp.py", "#/tmp/build_pysgpp/pysgpp.py", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/_pysgpp.so", "#/tmp/build_pysgpp/_pysgpp.so", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/pysgpp.py", "#/tmp/build_pysgpp/pysgpp.py", Copy("$TARGET", "$SOURCE"))

Help(opts.GenerateHelpText(env))
