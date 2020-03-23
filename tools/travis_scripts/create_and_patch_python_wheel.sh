#!/bin/bash

# exit when any command fails
set -e

# set paths
export LD_LIBRARY_PATH=$(pwd)/lib:$LD_LIBRARY_PATH

# create wheel and patch dependencies
python3 setup.py bdist_wheel
auditwheel repair dist/pysgpp-0.0.0-py3-none-any.whl --plat linux_x86_64

# install patched package
pip3 install wheelhouse/pysgpp-0.0.0-py3-none-linux_x86_64.whl

# test some examples
printf '\n\n\n~~~ Testing pip python package ~~~\n\n'
mkdir -p "$HOME/testing_pip_python_package"
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
  target_file=$HOME/testing_pip_python_package/${file_name}
  echo "... copying $file_name from" $source_file " into $target_file "
  cp ${source_file} ${target_file}
  echo "... switching into python test folder"
  cd "$HOME/testing_pip_python_package"
  echo "... running python3 ${file_name} "
  python3 ${file_name}
  echo "... switching back into SGpp root directory"
  cd -
done
