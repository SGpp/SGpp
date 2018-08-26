#!/bin/sh

scons -j $(nproc) VERBOSE=1 OPT=1 SG_ALL=0 SG_BASE=1 SG_OPTIMIZATION=1 SG_PYTHON=1 \
    USE_MATLAB=1 BUILD_STATICLIB=1 USE_PYTHON3_FOR_PYSGPP=1 RUN_PYTHON_TESTS=0 \
    CPPPATH=$REAL_HOME/local/share/matlab/R2017a/extern/include \
    LIBPATH=$REAL_HOME/local/share/matlab/R2017a/bin/glnxa64 \
    RPATH=$REAL_HOME/local/share/matlab/R2017a/bin/glnxa64
