#################################################################################
# Copyright (C) 2009-2011 Technische Universitaet Muenchen                      #
# This file is part of the SG++ project. For conditions of distribution and     #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp            #
#                                                                               #
# author Alexander Heinecke (Alexander.Heinecke@mytum.de)                       #
#################################################################################

###################################################################
# Needed Pathes
###################################################################
SRCDIR=./../../../src/sgpp
#only for extensions:
#####################
# Intel Array Building Blocks
ARBBINCLUDE = /opt/intel/arbb/1.0.0.030/include
ARBBLIB = /opt/intel/arbb/1.0.0.030/lib/intel64
# NVidia OpenCL
OCLINCLUDE = /opt/intel/intel-opencl-1.2-4.6.0.92/opencl-1.2-sdk-4.6.0.92/include/
OCLLIB = /opt/intel/intel-opencl-1.2-4.6.0.92/opencl-1.2-4.6.0.92/lib64/
#OCLINCLUDE = ${CUDATOOLKIT_HOME}/include
#OCLLIB = /opt/cray/nvidia/default/lib64
# Intel OpenCL
#IOCLINCLUDE = /usr/include
#IOCLLIB = /usr/lib64/OpenCL/vendors/intel
#IOCLINCLUDE = /etc/alternatives/opencl-headers
#IOCLLIB = /etc/alternatives/opencl-intel-runtime/lib64
IOCLINCLUDE = /opt/intel/intel-opencl-1.2-4.6.0.92/opencl-1.2-sdk-4.6.0.92/include/
IOCLLIB = /opt/intel/intel-opencl-1.2-4.6.0.92/opencl-1.2-4.6.0.92/lib64/
# AMD OpenCL
# AMDOCLINCLUDE = /opt/AMDAPP/include
# AMDOCLLIB = /opt/AMDAPP/lib/x86_64
AMDOCLINCLUDE = /lrz/sys/parallel/amdapp/2.7/include
AMDOCLLIB = /lrz/sys/parallel/amdapp/2.7/lib/x86_64
# Intel OpenCL, Windows
IOCLINCLUDEWIN = \"C:\Program Files (x86)\Intel\OpenCL SDK\3.0\include\"
IOCLLIBWIN = \"C:\Program Files (x86)\Intel\OpenCL SDK\3.0\lib\x64\OpenCL.lib\"

###################################################################
# Default Variables, overwirtten by CLI
###################################################################	
# use OpenMP Version 3
OMP=1
# use the TR1 Implementations for Hashmaps
TR1=0
# default compiler: g++; possible values: g++, icpc (Intel Compiler)
CC=g++
#CC=icpc
# vectorization option
#  sse3
#  sse4
#  avx
VEC=sse3
# extensions, manages extensions to be included, possible values (only when using Intel Compiler):
#       ArBB - Intel Array Building Blocks support
#       AMDOCLGPU - AMD GPU OpenCL support
#       INTELOCL - Intel CPU OpenCL support
#       INTELOCLGPU - Intel GPU OpenCL support
#       INTELOCLMIC - Intel MIC OpenCL support (Phi coprocessor)
#       NVOCL - NVIDIA OpenCL support
#       MIC_OFFLOAD - Intel MIC coprocessor in offload mode
#       MIC_NATIVE - compile for native execution on Intel MIC coprocessor
#       NO - no extensions, default
EXT=NVOCL
#	MPI - MPI support
MPI=0
# instances used to compile
JOBS=4
# Default residual threshold
SLE_RES_THRESH=-1.0
# Default number of parallel dimensions for the parallelization of the recursive up down scheme
UPDOWN_PARADIMS=4
# Compile for x86, but change vector width to match the one of MIC_NATIVE. This is needed for symmetric MPI execution. This option only has an effect for compilation with mpiicpc.
X86_MIC_SYMMETRIC=0
# Enable/Disable Operation Matrix Results
STORE_OP_MA=0
# Enable/Disable -ip -ipo
IPO=0


###################################################################
# Compiler Flags
###################################################################	
CFLAGS_GCC:=-Wno-deprecated-declarations -Wall -Werror -Wconversion -Wno-deprecated -Wno-long-long -pedantic -ansi -c -O3 -funroll-loops -fno-strict-aliasing -fPIC -mfpmath=sse -I$(SRCDIR)
LFLAGS_GCC:=-Wall -pedantic -ansi -O3

CFLAGS_ICC:=-Wall -Werror -wd1125 -Wconversion -Wno-deprecated -ipo -ip -ansi -ansi-alias -fp-speculation=safe -c -O3 -funroll-loops -fPIC -I$(SRCDIR)
LFLAGS_ICC:=-Wall -ipo -ip -ansi -O3

CFLAGS_ICL:=/Wall /Qipo /Qip /Oa /Qansi_alias /Qfp-speculation=safe /c /O3 /Qunroll-aggressive /I$(SRCDIR) /DUSETRONE /Qcxx-features /D_WIN32 /DNOMINMAX
LFLAGS_ICL:=/Wall /Qipo /Qip /Qansi_alias /O3

ifeq ($(IPO), 0)
CFLAGS_ICC:=-Wall -Werror -wd1125 -Wconversion -Wno-deprecated -ip -ansi -ansi-alias -fp-speculation=safe -c -O3 -funroll-loops -fPIC -I$(SRCDIR)
LFLAGS_ICC:=-Wall -ip -ansi -O3
endif

ifeq ($(EXT), MIC_NATIVE)
VEC:=""
endif

ifneq ($(EXT), MIC_OFFLOAD)
CFLAGS_ICC:=$(CFLAGS_ICC) -no-offload
LFLAGS_ICC:=$(LFLAGS_ICC) -no-offload
endif

ifeq ($(STORE_OP_MA), 1)
CFLAGS_ICC:=$(CFLAGS_ICC) -DSTORE_MATRIX=1
LFLAGS_ICC:=$(LFLAGS_ICC) -DSTORE_MATRIX=1

CFLAGS_GCC:=$(CFLAGS_GCC) -DSTORE_MATRIX=1
LFLAGS_GCC:=$(LFLAGS_GCC) -DSTORE_MATRIX=1
endif

ifeq ($(CC),g++)
CFLAGS:=$(CFLAGS_GCC)
LFLAGS:=$(LFLAGS_GCC)
#EXT=NO
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -fopenmp
LFLAGS:=$(LFLAGS) -fopenmp
endif
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -msse4.2
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) -mavx -D__USEAVX128__
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif
ifeq ($(VEC),avx2)
CFLAGS:=$(CFLAGS) -mavx2 -mfma
endif
ifeq ($(VEC),bd_avx128)
CFLAGS:=$(CFLAGS) -mavx -mfma4 -mxop -march=bdver1 -D__USEAVX128__
endif
ifeq ($(VEC),bd_avx)
CFLAGS:=$(CFLAGS) -mavx -mfma4 -mxop -march=bdver1
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -DUSEOCL_NVIDIA -fopenmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp -DUSEOCL_CPU
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLMIC)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -DUSEOCL_MIC
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) -I$(AMDOCLINCLUDE) -DUSEOCL -DUSEOCL_AMD -fopenmp
LFLAGS:=$(LFLAGS) -L$(AMDOCLLIB) -lOpenCL -fopenmp
endif
endif

ifeq ($(CC),icpc)
CFLAGS:=$(CFLAGS_ICC)
LFLAGS:=$(LFLAGS_ICC)
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -msse4.2
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) -mavx -D__USEAVX128__
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif
ifeq ($(VEC),avx2)
CFLAGS:=$(CFLAGS) -march=core-avx2 -fma
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -openmp
LFLAGS:=$(LFLAGS) -openmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -DUSEOCL_NVIDIA -openmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -openmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -openmp -DUSEOCL_CPU
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -openmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -openmp
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -openmp
endif
ifeq ($(EXT), INTELOCLMIC)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -DUSEOCL_MIC
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) -I$(AMDOCLINCLUDE) -DUSEOCL -DUSEOCL_AMD -openmp
LFLAGS:=$(LFLAGS) -L$(AMDOCLLIB) -lOpenCL -openmp
endif
ifeq ($(EXT), MIC_OFFLOAD)
CFLAGS:=$(CFLAGS) -no-ipo -no-ip -DUSEMIC -opt-report-phase:offload -openmp-report2 -offload-option,mic,compiler,\"-Wall -ansi -fPIC -c -ansi-alias -O3 -I$(SRCDIR) -openmp -fma -sox -mP2OPT_hlo_enable_all_mem_refs_prefetch -mP2OPT_hlo_pref_loop_based_prefetch=F\"
LFLAGS:=$(LFLAGS) -no-ipo -no-ip -DUSEMIC -offload-option,mic,link,\"--no-undefined\" 
XIAR_OPTS:=$(XIAR_OPTS) -qoffload-build
endif
ifeq ($(EXT), MIC_NATIVE)
CFLAGS:=$(CFLAGS) -mmic -DUSEMIC -fma -sox -mP2OPT_hlo_enable_all_mem_refs_prefetch -mP2OPT_hlo_pref_loop_based_prefetch=F -fp-speculation=fast
LFLAGS:=$(LFLAGS) -mmic
LIB_OPTS:=-mmic
endif
endif

ifeq ($(CC),icl)
CFLAGS:=$(CFLAGS_ICL)
LFLAGS:=$(LFLAGS_ICL)
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) /arch:SSE3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) /arch:SSE4.2
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) /arch:AVX /D__USEAVX128__
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) /arch:AVX
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) /Qopenmp
LFLAGS:=$(LFLAGS) /Qopenmp
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) /I$(OCLINCLUDE) /DUSEOCL /DUSEOCL_NVIDIA /Qopenmp
LFLAGS:=$(LFLAGS) /L$(OCLLIB) /Qpenmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) /I$(IOCLINCLUDEWIN) /DUSEOCL /DUSEOCL_INTEL /Qopenmp /DUSEOCL_CPU
LFLAGS:=$(LFLAGS) $(IOCLLIBWIN) /Qopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) /I$(IOCLINCLUDEWIN) /DUSEOCL /DUSEOCL_INTEL /Qopenmp
LFLAGS:=$(LFLAGS) $(IOCLLIBWIN) /Qopenmp
endif
ifeq ($(EXT), INTELOCLMIC)
CFLAGS:=$(CFLAGS) /I$(IOCLINCLUDE) /DUSEOCL /DUSEOCL_INTEL /DUSEOCL_MIC
LFLAGS:=$(LFLAGS) $(IOCLLIBWIN) /Qopenmp
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) /I$(IOCLINCLUDEWIN) /DUSEOCL /DUSEOCL_INTEL /Qopenmp
LFLAGS:=$(LFLAGS) /Qopenmp $(IOCLLIBWIN)
endif
endif

ifeq ($(CC),mpigxx)
CFLAGS:=$(CFLAGS_GCC)
LFLAGS:=$(LFLAGS_GCC)
CFLAGS:=$(CFLAGS) -DUSE_MPI
MPI=1
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -fopenmp
LFLAGS:=$(LFLAGS) -fopenmp
endif
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -msse4.2
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) -mavx -D__USEAVX128__
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif
ifeq ($(VEC),avx2)
CFLAGS:=$(CFLAGS) -mavx2 -fma
endif
ifeq ($(VEC),bd_avx128)
CFLAGS:=$(CFLAGS) -mavx -mfma4 -mxop -march=bdver1 -D__USEAVX128__
endif
ifeq ($(VEC),bd_avx)
CFLAGS:=$(CFLAGS) -mavx -mfma4 -mxop -march=bdver1
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -DUSEOCL_NVIDIA -fopenmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp -DUSEOCL_CPU
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLMIC)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -DUSEOCL_MIC
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) -I$(AMDOCLINCLUDE) -DUSEOCL -DUSEOCL_AMD -fopenmp
LFLAGS:=$(LFLAGS) -L$(AMDOCLLIB) -lOpenCL -fopenmp
endif
endif

ifeq ($(CC),CC)
CFLAGS:=-target=compute_node $(CFLAGS_GCC)
LFLAGS:=-target=compute_node $(LFLAGS_GCC)
CFLAGS:=$(CFLAGS) -DUSE_MPI
MPI=1
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -fopenmp
LFLAGS:=$(LFLAGS) -fopenmp
endif
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -msse4.2
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) -mavx -D__USEAVX128__
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif
ifeq ($(VEC),avx2)
CFLAGS:=$(CFLAGS) -mavx2 -fma
endif
ifeq ($(VEC),bd_avx128)
CFLAGS:=$(CFLAGS) -mavx -mfma4 -mxop -march=bdver1 -D__USEAVX128__
endif
ifeq ($(VEC),bd_avx)
CFLAGS:=$(CFLAGS) -mavx -mfma4 -mxop -march=bdver1
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -DUSEOCL_NVIDIA -fopenmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp -DUSEOCL_CPU
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) -I$(AMDOCLINCLUDE) -DUSEOCL -DUSEOCL_AMD -fopenmp
LFLAGS:=$(LFLAGS) -L$(AMDOCLLIB) -lOpenCL -fopenmp
endif
endif


ifeq ($(CC),mpiicpc)
CFLAGS:=$(CFLAGS_ICC)
LFLAGS:=$(LFLAGS_ICC)
CFLAGS:=$(CFLAGS) -DUSE_MPI
MPI=1
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -msse4.2
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif
ifeq ($(VEC),avx2)
CFLAGS:=$(CFLAGS) -march=core-avx2 -fma
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) -mavx -D__USEAVX128__
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -openmp
LFLAGS:=$(LFLAGS) -openmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -DUSEOCL_NVIDIA -fopenmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp -DUSEOCL_CPU
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -DUSEOCL_MIC
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) -I$(AMDOCLINCLUDE) -DUSEOCL -DUSEOCL_AMD -fopenmp
LFLAGS:=$(LFLAGS) -L$(AMDOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), MIC_OFFLOAD)
CFLAGS:=$(CFLAGS) -no-ipo -no-ip -DUSEMIC -opt-report-phase:offload -openmp-report2 -offload-option,mic,compiler,\"-Wall -ansi -fPIC -c -ansi-alias -O3 -I$(SRCDIR) -openmp -fma -sox -mP2OPT_hlo_enable_all_mem_refs_prefetch -mP2OPT_hlo_pref_loop_based_prefetch=F\"
LFLAGS:=$(LFLAGS) -no-ipo -no-ip -DUSEMIC -offload-option,mic,link,\"--no-undefined\" 
XIAR_OPTS:=$(XIAR_OPTS) -qoffload-build
endif
ifeq ($(EXT), MIC_NATIVE)
CFLAGS:=$(CFLAGS) -mmic -DUSEMIC -fma -sox -mP2OPT_hlo_enable_all_mem_refs_prefetch -mP2OPT_hlo_pref_loop_based_prefetch=F -fp-speculation=fast
LFLAGS:=$(LFLAGS) -mmic
LIB_OPTS:=-mmic
endif
ifeq ($(X86_MIC_SYMMETRIC), 1)
CFLAGS:=$(CFLAGS) -DX86_MIC_SYMMETRIC -DUSEMIC -DMIC_UNROLLING_WIDTH=96 -DMIC_UNROLLING_WIDTH_SP=192
endif
endif

ifeq ($(CC),mpiCC)
CFLAGS:=$(CFLAGS_ICC)
LFLAGS:=$(LFLAGS_ICC)
CFLAGS:=$(CFLAGS) -DUSE_MPI
MPI=1
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -msse4.2
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif
ifeq ($(VEC),avx2)
CFLAGS:=$(CFLAGS) -march=core-avx2 -fma
endif
ifeq ($(VEC),avx128)
CFLAGS:=$(CFLAGS) -mavx -D__USEAVX128__
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -openmp
LFLAGS:=$(LFLAGS) -openmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), NVOCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -DUSEOCL_NVIDIA -fopenmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCL)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp -DUSEOCL_CPU
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), INTELOCLGPU)
CFLAGS:=$(CFLAGS) -I$(IOCLINCLUDE) -DUSEOCL -DUSEOCL_INTEL -fopenmp
LFLAGS:=$(LFLAGS) -L$(IOCLLIB) -lOpenCL -fopenmp
endif
ifeq ($(EXT), AMDOCLGPU)
CFLAGS:=$(CFLAGS) -I$(AMDOCLINCLUDE) -DUSEOCL -DUSEOCL_AMD -fopenmp
LFLAGS:=$(LFLAGS) -L$(AMDOCLLIB) -lOpenCL -fopenmp
endif
endif

CFLAGS:=$(CFLAGS) -DDEFAULT_RES_THRESHOLD=$(SLE_RES_THRESH)
CFLAGS:=$(CFLAGS) -DTASKS_PARALLEL_UPDOWN=$(UPDOWN_PARADIMS)
CFLAGS:=$(CFLAGS) -DSG_PARALLEL

###################################################################
# Builds a lib containing all SG Algorithms
###################################################################	
default:
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/sgpplib_gcc
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc" "EXT=$(EXT)"
endif
ifeq ($(CC),opencc)
	mkdir -p tmp/build_native/sgpplib_opencc
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_opencc "CC=/opt/x86_open64-4.2.5.2/bin/$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_opencc" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/sgpplib_icc
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc" "EXT=$(EXT)" "XIAR_OPTS=$(XIAR_OPTS)" "LIB_OPTS=$(LIB_OPTS)"
endif
ifeq ($(CC),icl)
	mkdir -p tmp/build_native/sgpplib_icl
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_icl "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icl" "EXT=$(EXT)"
endif
ifeq ($(CC),CC)
	mkdir -p tmp/build_native/sgpplib_cray
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_cray "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_cray" "EXT=$(EXT)" "MPI=$(MPI)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/sgpplib_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc" "EXT=$(EXT)" "MPI=$(MPI)" "XIAR_OPTS=$(XIAR_OPTS)" "LIB_OPTS=$(LIB_OPTS)"
endif
ifeq ($(CC),mpigxx)
	mkdir -p tmp/build_native/sgpplib_mpigxx
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_mpigxx "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpigxx" "EXT=$(EXT)" "MPI=$(MPI)"
endif
ifeq ($(CC),mpiCC)
	mkdir -p tmp/build_native/sgpplib_ibm
	make -j $(JOBS) -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_ibm "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_ibm" "EXT=$(EXT)" "MPI=$(MPI)"
endif


###################################################################
# Builds a Balck Scholes Solver
###################################################################	
BSSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSSolver_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolver --directory=./tmp/build_native/BSSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSSolver_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolver --directory=./tmp/build_native/BSSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSSolver_ICC" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/BSSolver_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolverMPI --directory=./tmp/build_native/BSSolver_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=BSSolver_ICC_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),mpigxx)
	mkdir -p tmp/build_native/BSSolver_mpigxx
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolverMPI --directory=./tmp/build_native/BSSolver_mpigxx "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpigxx.a" "BINNAME=BSSolver_GCC_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),CC)
	mkdir -p tmp/build_native/BSSolver_cray
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolverMPI --directory=./tmp/build_native/BSSolver_cray "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_cray.a" "BINNAME=BSSolver_CRAY_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiCC)
	mkdir -p tmp/build_native/BSSolver_ibm
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolverMPI --directory=./tmp/build_native/BSSolver_ibm "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_ibm.a" "BINNAME=BSSolver_IBM_MPI" "EXT=$(EXT)"
endif

###################################################################
# Builds a Black Scholes Solver with Stretching
###################################################################	
BSSolverWithStretching: default

ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSSolverWithStretching_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolverWithStretching --directory=./tmp/build_native/BSSolverWithStretching_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSSolverWithStretching_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSSolverWithStretching_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeBlackScholesSolverWithStretching --directory=./tmp/build_native/BSSolverWithStretching_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSSolverWithStretching_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Heston Solver
####################################################################  

HestonSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HestonSolver_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeHestonSolver --directory=./tmp/build_native/HestonSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HestonSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HestonSolver_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeHestonSolver --directory=./tmp/build_native/HestonSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HestonSolver_ICC" "EXT=$(EXT)"
endif
#ifeq ($(CC),mpiicpc)   
#   Not implemented     
#endif

###################################################################
# Builds a Hull White Solver
###################################################################	
HWSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HWSolver_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeHullWhiteSolver --directory=./tmp/build_native/HWSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HWSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HWSolver_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeHullWhiteSolver --directory=./tmp/build_native/HWSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HWSolver_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Hull White combine Black Scholes Solver
###################################################################	
BSHWSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSHWSolver_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeBSHWSolver --directory=./tmp/build_native/BSHWSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSHWSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSHWSolver_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeBSHWSolver --directory=./tmp/build_native/BSHWSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSHWSolver_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple Heat Equation Solver
###################################################################	
HESolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HESolver_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolver --directory=./tmp/build_native/HESolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HESolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HESolver_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolver --directory=./tmp/build_native/HESolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HESolver_ICC" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/HESolver_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolverMPI --directory=./tmp/build_native/HESolver_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=HESolver_ICC_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),mpigxx)
	mkdir -p tmp/build_native/HESolver_mpigxx
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolverMPI --directory=./tmp/build_native/HESolver_mpigxx "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpigxx.a" "BINNAME=HESolver_GCC_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiCC)
	mkdir -p tmp/build_native/HESolver_ibm
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolverMPI --directory=./tmp/build_native/HESolver_ibm "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_ibm.a" "BINNAME=HESolver_IBM_MPI" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple Heat Equation Solver (rotating Laser test case)
###################################################################	
LaserHESolver2D: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/LaserHESolver2D_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeLaserHeatEquationSolver --directory=./tmp/build_native/LaserHESolver2D_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=LaserHESolver2D_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/LaserHESolver2D_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeLaserHeatEquationSolver --directory=./tmp/build_native/LaserHESolver2D_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=LaserHESolver2D_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple Heat Equation Solver with Stretching
###################################################################	
HESolverWithStretching: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HESolverWithStretching_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolverWithStretching --directory=./tmp/build_native/HESolverWithStretching_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HESolverWithStretching_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HESolverWithStretching_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeHeatEquationSolverWithStretching --directory=./tmp/build_native/HESolverWithStretching_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HESolverWithStretching_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a ClassifyBenchmark Application
###################################################################	
ClassifyBenchmark: default
ifeq ($(CC),CC)
	mkdir -p tmp/build_native/ClassifyBenchmark_cray
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmarkMPI --directory=./tmp/build_native/ClassifyBenchmark_cray "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_cray.a" "BINNAME=ClassifyBenchmark_CRAY_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/ClassifyBenchmark_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=ClassifyBenchmark_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),opencc)
	mkdir -p tmp/build_native/ClassifyBenchmark_opencc
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_opencc "CC=/opt/x86_open64-4.2.5.2/bin/$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_opencc.a" "BINNAME=ClassifyBenchmark_OPENCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/ClassifyBenchmark_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=ClassifyBenchmark_ICC" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/ClassifyBenchmark_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmarkMPI --directory=./tmp/build_native/ClassifyBenchmark_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=ClassifyBenchmark_ICC_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),mpigxx)
	mkdir -p tmp/build_native/ClassifyBenchmark_mpigxx
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmarkMPI --directory=./tmp/build_native/ClassifyBenchmark_mpigxx "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpigxx.a" "BINNAME=ClassifyBenchmark_GCC_MPI" "EXT=$(EXT)"
endif
ifeq ($(CC),icl)
	mkdir -p tmp/build_native/ClassifyBenchmark_icl
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_icl "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icl.lib" "BINNAME=ClassifyBenchmark_ICL.exe" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiCC)
	mkdir -p tmp/build_native/ClassifyBenchmark_ibm
	make -j $(JOBS) -f ./../../../src/makefileNativeClassifyBenchmarkMPI --directory=./tmp/build_native/ClassifyBenchmark_ibm "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_ibm.a" "BINNAME=ClassifyBenchmark_IBM_MPI" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple LaplaceTest App
###################################################################	
LaplaceTest: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/LaplaceTest_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeLaplaceTest --directory=./tmp/build_native/LaplaceTest_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=LaplaceTest_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/LaplaceTest_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeLaplaceTest --directory=./tmp/build_native/LaplaceTest_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=LaplaceTest_ICC" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/LaplaceTest_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileNativeLaplaceTestMPI --directory=./tmp/build_native/LaplaceTest_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=LaplaceTest_MPIICC" "EXT=$(EXT)"
endif


###################################################################
# Builds a simple LTwoDotProductTest App
###################################################################
LTwoDotProductTest: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/LTwoDotProductTest_gcc
	make -j $(JOBS) -f ./../../../src/makefileNativeLTwoDotProductTest --directory=./tmp/build_native/LTwoDotProductTest_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=LTwoDotProductTest_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/LTwoDotProductTest_icc
	make -j $(JOBS) -f ./../../../src/makefileNativeLTwoDotProductTest --directory=./tmp/build_native/LTwoDotProductTest_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=LTwoDotProductTest_ICC" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/LTwoDotProductTest_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileNativeLTwoDotProductTestMPI --directory=./tmp/build_native/LTwoDotProductTest_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=LTwoDotProductTest_MPIICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple LaplaceLTwoDotProductTest App
###################################################################
LaplaceLTwoDotProductTest: default
#ifeq ($(CC),g++)
#	mkdir -p tmp/build_native/LaplaceLTwoDotProductTest_gcc
#	make -j $(JOBS) -f ./../../../src/makefileNativeLaplaceLTwoDotProductTest --directory=./tmp/build_native/LaplaceLTwoDotProductTest_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=LaplaceLTwoDotProductTest_GCC" "EXT=$(EXT)"
#endif
#ifeq ($(CC),icpc)
#	mkdir -p tmp/build_native/LaplaceLTwoDotProductTest_icc
#	make -j $(JOBS) -f ./../../../src/makefileNativeLaplaceLTwoDotProductTest --directory=./tmp/build_native/LaplaceLTwoDotProductTest_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=LaplaceLTwoDotProductTest_ICC" "EXT=$(EXT)"
#endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/LaplaceLTwoDotProductTest_mpiicc
	make -j $(JOBS) -f ./../../../src/makefileNativeLaplaceLTwoDotProductTestMPI --directory=./tmp/build_native/LaplaceLTwoDotProductTest_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=LaplaceLTwoDotProductTest_MPIICC" "EXT=$(EXT)"
endif


###################################################################
# test Balck Scholes Solver
###################################################################	
		
test_BS_1d:
	cd bin; \
	./copyBSSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolver/1d; \
	./test_BSSolver_1d.sh;
	
test_BS_2d:
	cd bin; \
	./copyBSSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolver/2d; \
	./test_BSSolver_2d.sh;
		
test_BS_3d:
	cd bin; \
	./copyBSSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolver/3d; \
	./test_BSSolver_3d.sh;
	
test_BS_all: test_BS_1d test_BS_2d test_BS_3d
	echo "executed all BS tests!"

###################################################################
# test Black Scholes Solver with Stretching
###################################################################	
		
test_BSS_1d:
	cd bin; \
	./copyBSSolverWithStretchingToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolverWithStretching/1d; \
	./test_BSSolverWithStretching_1d.sh;
	
test_BSS_2d:
	cd bin; \
	./copyBSSolverWithStretchingToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolverWithStretching/2d; \
	./test_BSSolverWithStretching_2d.sh;
		
test_BSS_3d:
	cd bin; \
	./copyBSSolverWithStretchingToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolverWithStretching/3d; \
	./test_BSSolverWithStretching_3d.sh;
	
test_BSS_all: test_BSS_1d test_BSS_2d test_BSS_3d
	echo "executed all BS tests!"
	
###################################################################
# test Heston Solver    
# ###################################################################                     

test_Heston_1d:
	cd bin; \
	./copyHestonSolverToTest.sh; \
	cd ./../tests/CPP_Apps/HestonSolver/1d; \
	./test_HestonSolver_1d.sh;

test_Heston_2d:
	cd bin; \
	./copyHestonSolverToTest.sh; \
	cd ./../tests/CPP_Apps/HestonSolver/2d; \
	./test_HestonSolver_2d.sh;
                                                                        
test_Heston_all: test_Heston_1d test_Heston_2d
	echo "executed all Heston tests!"
                                                                                 	
###################################################################
# test Combined Hull Wihte Solver Solver
###################################################################			

test_BSHW:
	cd bin; \
	./copyBSHWSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSHWSolver; \
	./test_BSHWSolver_2d_cart.sh;
		
###################################################################
# test ClassifyBenchmark
###################################################################	

test_ClassifyBenchmark:
	cd bin; \
	./copyClassifyBenchmarkToTest.sh; \
	cd ./../tests/CPP_Apps/ClassifyBenchmark; \
	./test_ClassifyBenchmark.sh;

clean:
	rm -rfv tmp/build_native
