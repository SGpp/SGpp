#!/bin/sh

cd bin
./copyBSSolverToTest.sh
cd ..
cd ./tests/BlackScholes_TestData
./test_BSSolver.sh
cd ./../..
