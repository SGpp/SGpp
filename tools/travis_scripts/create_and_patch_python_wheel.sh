#!/bin/bash

# exit when any command fails
set -e

pushd ../..

# set paths
export LD_LIBRARY_PATH=$(pwd)/lib/sgpp:$LD_LIBRARY_PATH

# create wheel and patch dependencies
python3 setup.py bdist_wheel
auditwheel repair dist/pysgpp-3.3.0-py3-none-any.whl --plat linux_x86_64

# install patched package
pip3 install wheelhouse/pysgpp-3.3.0-py3-none-linux_x86_64.whl

# test base example
python3 base/examples/quadrature.py