#############################################################################
# This file is part of sgpp, a program package making use of spatially      #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

SRCDIR=./../../../src/sgpp
###################################################################
# Default Variables, overwirtten by CLI
###################################################################	
# use OpenMP Version 3
OMP=0
# use the TR1 Implementations for Hashmaps
TR1=0
# default compiler: g++; possible values: g++, icpc (Intel Compiler)
CC=g++
# vectorization option for intel compiler
VEC=sse3

###################################################################
# Compiler Flags
###################################################################	
CFLAGS_GCC:=-Wall -pedantic -ansi -c -Wno-long-long -fno-strict-aliasing -O3 -funroll-loops -ffloat-store -I$(SRCDIR)
LFLAGS_GCC:=-Wall -pedantic -ansi -O3

CFLAGS_ICC:=-Wall -ansi -c -fno-strict-aliasing -ipo -ip -ansi-alias -O3 -funroll-loops -I$(SRCDIR) -DUSEICCINTRINSICS
LFLAGS_ICC:=-Wall -ansi -O3 -static-intel -ipo -ip

ifeq ($(CC),g++)
CFLAGS:=$(CFLAGS_GCC)
LFLAGS:=$(LFLAGS_GCC)
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -fopenmp
LFLAGS:=$(LFLAGS) -fopenmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
endif
ifeq ($(CC),icpc)
CFLAGS:=$(CFLAGS_ICC)
LFLAGS:=$(LFLAGS_ICC)
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -xSSE3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -xSSE4.2
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -xAVX
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -openmp
LFLAGS:=$(LFLAGS) -openmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
endif

###################################################################
# Builds a lib containing all SG Algorithms
###################################################################	
default:
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/sgpplib_gcc
	make -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/sgpplib_icc
	make -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a"
endif

###################################################################
# Builds a Balck Scholes Solver
###################################################################	
BSSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSSolver_gcc
	make -f ./../../../src/makefileNativeBlackScholesSolver --directory=./tmp/build_native/BSSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSSolver_GCC"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSSolver_icc
	make -f ./../../../src/makefileNativeBlackScholesSolver --directory=./tmp/build_native/BSSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSSolver_ICC"
endif

###################################################################
# Builds a Hull White Solver
###################################################################	
HWSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HWSolver_gcc
	make -f ./../../../src/makefileNativeHullWhiteSolver --directory=./tmp/build_native/HWSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HWSolver_GCC"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HWSolver_icc
	make -f ./../../../src/makefileNativeHullWhiteSolver --directory=./tmp/build_native/HWSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HWSolver_ICC"
endif

###################################################################
# Builds a Hull White combine Black Scholes Solver
###################################################################	
BSHWSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSHWSolver_gcc
	make -f ./../../../src/makefileNativeBSHWSolver --directory=./tmp/build_native/BSHWSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSHWSolver_GCC"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSHWSolver_icc
	make -f ./../../../src/makefileNativeBSHWSolver --directory=./tmp/build_native/BSHWSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSHWSolver_ICC"
endif

###################################################################
# Builds a simple Heat Equation Solver
###################################################################	
HESolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HESolver_gcc
	make -f ./../../../src/makefileNativeHeatEquationSolver --directory=./tmp/build_native/HESolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HESolver_GCC"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HESolver_icc
	make -f ./../../../src/makefileNativeHeatEquationSolver --directory=./tmp/build_native/HESolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HESolver_ICC"
endif

###################################################################
# Builds a ClassifyBenchmark Application
###################################################################	
ClassifyBenchmark: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/ClassifyBenchmark_gcc
	make -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=ClassifyBenchmark_GCC"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/ClassifyBenchmark_icc
	make -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=ClassifyBenchmark_ICC"
endif

###################################################################
# Builds a Up/Down Test Application
###################################################################	
UpDownTest: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/UpDownTest_gcc
	make -f ./../../../src/makefileUpDownTest --directory=./tmp/build_native/UpDownTest_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=UpDownTest_GCC"
endif

###################################################################
# Builds a Refine/Coarsen Test Application
###################################################################	
RefineCoarsenTest: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/RefineCoarsen_gcc
	make -f ./../../../src/makefileRefineCoarsenTest --directory=./tmp/build_native/RefineCoarsen_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=RefineCoarsen_gcc"
endif
		
clean:
	rm -rdfv tmp/build_native
