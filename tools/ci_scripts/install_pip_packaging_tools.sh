#!/bin/bash
set -e

export PIP_BREAK_SYSTEM_PACKAGES=1
# upgrade pip
sudo -H pip3 install --upgrade pip
# install setuptools/wheel
sudo -H pip3 install setuptools wheel scipy numpy jinja2

# install patchelf
git clone https://github.com/NixOS/patchelf.git
pushd patchelf
git checkout tags/0.14.5
./bootstrap.sh
./configure
make
sudo make install
popd

# install auditwheel
sudo -H pip3 install auditwheel
