#!/bin/sh
module load boost-1.69.0-gcc-4.8.5-gbslxeg

scons -j $(nproc) VERBOSE=1 OPT=1 SG_ALL=0 SG_BASE=1 SG_OPTIMIZATION=1 SG_PYTHON=0 \
    USE_MATLAB=1 BUILD_STATICLIB=1 USE_PYTHON3_FOR_PYSGPP=0 RUN_PYTHON_TESTS=0 \
    CPPPATH=/home/portfolio/public/MATLAB/R2017b/extern/include \
    LIBPATH=/home/portfolio/public/MATLAB/R2017b/bin/glnxa64 \
    RPATH=/home/portfolio/public/MATLAB/R2017b/bin/glnxa64
