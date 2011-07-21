#!/bin/sh

cd cart
./test_BSSolver_3d_cart.sh
cd ..
cd log
./test_BSSolver_3d_log.sh
./test_BSSolver_3d_PAT.sh
cd ..
