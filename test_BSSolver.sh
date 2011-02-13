#!/bin/sh

cd bin
./copyBSSolverToTest.sh
cd ..
cd ./tests/BlackScholes_TestData
./test_BSHWSolver.sh
./test_BSSolver.sh
cd ./../..
