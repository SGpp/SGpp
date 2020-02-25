#!/bin/bash

# exit when any command fails
set -e

printf '\n\n\n~~~ Creating debian packages ~~~\n\n'
cd tools/create_release/deb_package
./create_deb.sh testing

printf '\n\n\n~~~ Installing debian packages ~~~\n\n'
sudo dpkg -i libsgpp-test-package_0.0-0.deb

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
  echo "... running ./${exec_name} "
  ./${exec_name}
  echo "... switching back into SGpp root directory"
  cd -
done

printf '\n\n\n~~~ Testing debian python package ~~~\n\n'
sudo dpkg -i libsgpp-python-test-package_0.0-0.deb
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

printf '\n\n\n~~~ Testing debian java package ~~~\n\n'
sudo dpkg -i libsgpp-java-test-package_0.0-0.deb
mkdir -p "$HOME/testing_java_package"
declare -a java_examples_to_test=(
                                    "base/examples/tutorial.java"
                                    "base/examples/refinement.java"
                                )
arraylength=${#java_examples_to_test[@]}
for (( i=1; i<${arraylength}+1; i++ ));
do
  echo ""
  echo "-------------------------------------------------------------"
  echo "Java examples test" $i "of" ${arraylength} "(to be copied from " ${java_examples_to_test[$i-1]} "into test folder) : "
  source_file=${java_examples_to_test[$i-1]}
  file_name=$(basename $source_file)
  target_file=$HOME/testing_java_package/${file_name}
  echo "... copying $file_name from" $source_file " into $target_file "
  cp ${source_file} ${target_file}
  echo "... switching into java test folder"
  cd "$HOME/testing_java_package"
  exec_name=${file_name%.*}
  echo "... building ${target_file} into $exec_name "
  export CLASSPATH=/usr/share/java/jsgpp.jar:$(pwd)
  javac ${file_name}
  echo "... running java ${exec_name} "
  java ${exec_name}
  echo "... switching back into SGpp root directory"
  cd -
done
echo ""
echo "-------------------------------------------------------------"
echo "Last java examples test (needs two files from optimzation/examples) : "
echo "... copying optimization.java from base/examples/optimization.java into $HOME/testing_java_package/optimization.java"
cp base/examples/optimization.java $HOME/testing_java_package/optimization.java
echo "... copying ExampleFunction.java from base/examples/ExampleFunction.java into $HOME/testing_java_package/ExampleFunction.java"
cp base/examples/ExampleFunction.java $HOME/testing_java_package/ExampleFunction.java
echo "... switching into java test folder"
cd "$HOME/testing_java_package"
echo "... building optimization.java and ExampleFunction.java into optimization"
export CLASSPATH=/usr/share/java/jsgpp.jar:$(pwd)
javac optimization.java ExampleFunction.java
echo "... running java optimization "
java optimization
echo "... switching back into SGpp root directory"
cd -


# Cleanup
sudo dpkg --remove libsgpp libsgpp-python libsgpp-java
