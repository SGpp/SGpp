#!/bin/sh
module load boost-1.69.0-gcc-4.8.5-gbslxeg
module load python-2.7.15-gcc-4.8.5-7ok6o2a
module load scons-2.5.1-gcc-4.8.5-doh7q4i
module load comp/gcc/8.2.0

scons -j $(nproc) VERBOSE=1 OPT=1 SG_ALL=0 SG_BASE=1 SG_OPTIMIZATION=1 SG_PYTHON=0 \
    USE_MATLAB=1 BUILD_STATICLIB=1 USE_PYTHON3_FOR_PYSGPP=0 RUN_PYTHON_TESTS=0 \
    CPPPATH=/home/ww_matlab/MATLAB/R2019b/extern/include \
    LIBPATH=/home/ww_matlab/MATLAB/R2019b/bin/glnxa64 \
    RPATH=/home/ww_matlab/MATLAB/R2019b/bin/glnxa64
