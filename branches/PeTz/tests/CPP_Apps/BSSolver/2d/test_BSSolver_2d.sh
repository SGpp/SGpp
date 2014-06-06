#!/bin/sh

cd cart
./test_BSSolver_2d_cart.sh
cd ..
cd log
./test_BSSolver_2d_log.sh
./test_BSSolver_2d_PAT.sh
cd ..
