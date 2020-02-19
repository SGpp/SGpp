#!/bin/bash

# exit when any command fails
set -e

printf '\n~~~ building sg++ ~~~\n\n'
scons -j4 SG_JAVA=0 RUN_BOOST_TESTS=0 RUN_PYTHON_TESTS=0 CHECK_STYLE=0
cd tools/create_release/deb_package

# Create dependency lists 
# IMPORTANT: list needs to be comma delimited!
echo -e 'zlib1g, build-essential' > deps_sgpp.tmp
echo -e '' > deps_pysgpp.tmp

# Build debian package
printf '\n~~~ Building Debian package ~~~\n\n'
python3 create_deb.py deps_sgpp.tmp deps_pysgpp.tmp automated


printf '\n~~~ Installing debian package ~~~\n\n'
sudo dpkg -i libsgpp-test-package_0.0-0.deb
sudo dpkg -i libsgpp-python-test-package_0.0-0.deb

printf '\n~~~ Testing debian cpp package ~~~\n\n'
cd ../../..
mkdir -p "$HOME/testing_package"
cp "base/examples/quadrature.cpp" "$HOME/testing_package/quadrature.cpp"
cp "base/examples/quadrature.py" "$HOME/testing_package/quadrature.py"
cd "$HOME/testing_package"
g++ quadrature.cpp -std=c++11 -o quad -l sgppbase
./quad > cpp_output.txt
cat cpp_output

printf '\n~~~ Testing debian python package ~~~\n\n'
python3 quadrature.py

# Clean-up (for debugging, the lines below should be commented out)
# rm deps_sgpp.tmp
# rm deps_pysgpp.tmp
