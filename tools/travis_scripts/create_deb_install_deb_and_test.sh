#!/bin/bash

# exit when any command fails
set -e

printf '\n\n\n~~~ Creating debian packages ~~~\n\n'
cd tools/create_release/deb_package
./create_deb.sh testing

printf '\n\n\n~~~ Installing debian packages ~~~\n\n'
sudo dpkg -i libsgpp-test-package_0.0-0.deb
sudo dpkg -i libsgpp-python-test-package_0.0-0.deb

printf '\n\n\n~~~ Testing debian cpp package ~~~\n\n'
cd ../../..
mkdir -p "$HOME/testing_package"
declare -a cpp_examples_to_test=(
                                    "base/examples/quadrature.cpp"
                                    "base/examples/tutorial.cpp"
                                    "base/examples/gridTExample.cpp"
                                    "base/examples/JSONExample.cpp"
                                    "base/examples/gridInteractionExample.cpp"
                                    "base/examples/unserializeGrid.cpp"
                                    "base/examples/gridTypes.cpp"
                                    "base/examples/dataVectorSerializeDemo.cpp"
                                    "base/examples/refinement.cpp"
                                    "base/examples/predictiveRefinement.cpp"
                                    "solver/examples/fistaExample.cpp"
                                    "optimization/examples/optimization.cpp"
                                    "optimization/examples/constrainedOptimization.cpp"
                                    "combigrid/examples/combigrid.cpp" 
                                )
arraylength=${#cpp_examples_to_test[@]}
for (( i=1; i<${arraylength}+1; i++ ));
do
  echo ""
  echo "-------------------------------------------------------------"
  echo "CPP examples test" $i "of" ${arraylength} "(to be copied from " ${cpp_examples_to_test[$i-1]} "into test folder) : "
  source_file=${cpp_examples_to_test[$i-1]}
  file_name=$(basename $source_file)
  target_file=$HOME/testing_package/${file_name}
  echo "... copying $file_name from" $source_file " into $target_file "
  cp ${source_file} ${target_file}
  exec_name=${file_name%.*}
  echo "... switching into test folder"
  cd "$HOME/testing_package"
  echo "... building ${target_file} into $exec_name "
  g++ ${file_name} -std=c++11 -o ${exec_name} -l sgppbase -l sgppsolver -l sgppoptimization -lsgppcombigrid
  echo "... running ./$exec_name "
  ./${exec_name}
  echo "... switching back into SGpp root directory"
  cd -
done

printf '\n\n\n~~~ Testing debian python package ~~~\n\n'
mkdir -p "$HOME/testing_python_package"
declare -a python_examples_to_test=(
                                    "base/examples/quadrature.py"
                                    "base/examples/tutorial.py"
                                    "base/examples/gridTExample.py"
                                    "base/examples/dataVectorSerializeDemo.py"
                                    "base/examples/refinement.py"
                                    "base/examples/predictiveRefinement.py"
                                    "base/examples/predictiveANOVARefinement.py"
                                    "base/examples/subspaceRefinement.py"
                                    "optimization/examples/optimization.py"
                                    "optimization/examples/fuzzy.py"
                                    "combigrid/examples/combigrid.py"
                                    "pde/examples/LTwoDotTest.py"
                                )
arraylength=${#python_examples_to_test[@]}
for (( i=1; i<${arraylength}+1; i++ ));
do
  echo ""
  echo "-------------------------------------------------------------"
  echo "Python examples test" $i "of" ${arraylength} "(to be copied from " ${python_examples_to_test[$i-1]} "into test folder) : "
  source_file=${python_examples_to_test[$i-1]}
  file_name=$(basename $source_file)
  target_file=$HOME/testing_python_package/${file_name}
  echo "... copying $file_name from" $source_file " into $target_file "
  cp ${source_file} ${target_file}
  echo "... switching into python test folder"
  cd "$HOME/testing_python_package"
  echo "... running python3 ${file_name} "
  python3 ${file_name}
  echo "... switching back into SGpp root directory"
  cd -
done

# Cleanup
sudo dpkg -r libsgpp-test-package_0.0-0.deb
sudo dpkg -r libsgpp-python-test-package_0.0-0.deb
