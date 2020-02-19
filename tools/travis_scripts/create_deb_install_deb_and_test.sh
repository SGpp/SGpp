#!/bin/bash
printf '\nWARNING: This is only meant to run on travis! If you are running this on your own you probably want to run create_deb.sh instead!\n'
sleep 7
printf '\nOkay you are on your own now!\n'
sleep 1

# Build SG++ in the background
printf '\n~~~ building sg++ ~~~\n\n'
scons -j4 SG_JAVA=0 RUN_BOOST_TESTS=0 RUN_PYTHON_TESTS=0
cd tools/create_release/deb_package

# Create dependency lists 
# IMPORTANT: list needs to be comma delimited!
echo -e 'zlib1g, build-essential' > deps_sgpp.tmp
echo -e '' > deps_pysgpp.tmp

# Build debian package
printf '\n~~~ Building Debian package ~~~\n\n'
python3 create_deb.py deps_sgpp.tmp deps_pysgpp.tmp automated


printf '\n~~~ Installing debian package ~~~\n\n'
dpkg -i libsgpp-test-package_0.0-0.deb
dpkg -i libsgpp-python-test-package_0.0-0.deb

# Clean-up (for debugging, the lines below should be commented out)
# rm deps_sgpp.tmp
# rm deps_pysgpp.tmp
