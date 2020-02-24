#!/bin/bash
set -e

# upgrade pip
sudo -H pip3 install --upgrade pip
# install setuptools/wheel
sudo -H pip3 install setuptools wheel scipy numpy jinja2

# install patchelf
git clone https://github.com/NixOS/patchelf.git
pushd patchelf
git checkout tags/0.10
./bootstrap.sh
./configure
make
sudo make install
popd

# install auditwheel
sudo -H pip3 install auditwheel
