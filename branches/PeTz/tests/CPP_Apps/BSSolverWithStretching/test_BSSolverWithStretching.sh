#!/bin/sh

cd 1d
./test_BSSolverWithStretching_1d.sh
cd ..
cd 2d
./test_BSSolverWithStretching_2d.sh
cd ..
cd 3d
./test_BSSolverWithStretching_3d.sh
cd ..
