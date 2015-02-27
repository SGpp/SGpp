#!/bin/sh

cd cart
./test_BSSolver_1d_cart.sh
cd ..
cd log
./test_BSSolver_1d_log.sh
./test_BSSolver_1d_PAT.sh
cd ..
