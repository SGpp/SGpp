#!/bin/bash

# Build SG++ in the background
printf '\n~~~ Building SG++ ~~~\n\n'
cd ../../..
scons -j4 SG_JAVA=0 RUN_BOOST_TESTS=0 CHECK_STYLE=0
cd -

# Create dependency lists 
# IMPORTANT: list needs to be comma delimited!
echo -e 'zlib1g, build-essential' > deps_sgpp.tmp
echo -e '' > deps_pysgpp.tmp

# Build debian package
printf '\n~~~ Building Debian package ~~~\n\n'
python3 create_deb.py deps_sgpp.tmp deps_pysgpp.tmp

# Clean-up (for debugging, the lines below should be commented out)
# rm deps_sgpp.tmp
# rm deps_pysgpp.tmp
