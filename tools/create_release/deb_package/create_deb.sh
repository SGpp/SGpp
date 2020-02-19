#!/bin/bash
set -e

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
# If $1 is empty, create_deb.py will query the user for version numbers, author names, etc
# Use that for the actual release since you want to add the version numbers
# Otherwise, use some text like "./create_deb.sh testing" for $1 if you just want a package for automated testing
python3 create_deb.py deps_sgpp.tmp deps_pysgpp.tmp $1

# Clean-up (for debugging, the lines below should be commented out)
# rm deps_sgpp.tmp
# rm deps_pysgpp.tmp
