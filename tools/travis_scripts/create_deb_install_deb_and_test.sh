#!/bin/bash

# exit when any command fails
set -e

cd tools/create_release/deb_package
./create_deb.sh testing

printf '\n~~~ Installing debian package ~~~\n\n'
sudo dpkg -i libsgpp-test-package_0.0-0.deb
sudo dpkg -i libsgpp-python-test-package_0.0-0.deb

printf '\n~~~ Testing debian cpp package ~~~\n\n'
cd ../../..
mkdir -p "$HOME/testing_package"
cp "base/examples/quadrature.cpp" "$HOME/testing_package/quadrature.cpp"
cp "base/examples/quadrature.py" "$HOME/testing_package/quadrature.py"
cd "$HOME/testing_package"
printf "\nBuilding quadrature cpp example in directory $(pwd)"
g++ quadrature.cpp -std=c++11 -o quad -l sgppbase
./quad 

printf '\n~~~ Testing debian python package ~~~\n\n'
printf "\nBuilding quadrature cpp example in directory $(pwd)"
python3 quadrature.py

# Clean-up (for debugging, the lines below should be commented out)
# rm deps_sgpp.tmp
# rm deps_pysgpp.tmp
