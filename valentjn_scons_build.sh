#!/bin/sh

scons -j $(nproc) VERBOSE=1 OPT=1 SG_ALL=0 SG_BASE=1 SG_OPTIMIZATION=1 SG_PYTHON=1 \
    USE_EIGEN=1 USE_MATLAB=1 BUILD_STATICLIB=1 \
    CPPPATH=/home/valentjn/local/share/matlab/R2016a/extern/include \
    LIBPATH=/home/valentjn/local/share/matlab/R2016a/bin/glnxa64 \
    RPATH=/home/valentjn/local/share/matlab/R2016a/bin/glnxa64
