#!/bin/sh
module load python/2.7.15
module load gcc/5.2.0
module load boost/1.52.0
module list

scons -j $(nproc) VERBOSE=1 OPT=1 SG_ALL=0 SG_BASE=1 SG_OPTIMIZATION=1 SG_PYTHON=0 \
    USE_MATLAB=1 BUILD_STATICLIB=1 USE_PYTHON3_FOR_PYSGPP=0 RUN_PYTHON_TESTS=0 \
    CPPPATH=/home/portfolio/schober/opt/MATLAB/R2017b/extern/include \
    LIBPATH=/home/portfolio/schober/opt/MATLAB/R2017b/bin/glnxa64 \
    RPATH=/home/portfolio/schober/opt/MATLAB/R2017b/bin/glnxa64
