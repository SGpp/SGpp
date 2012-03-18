#!/bin/sh

cd cart
./test_HestonSolver_1d_cart.sh
cd ..
cd log
./test_HestonSolver_1d_log.sh
cd ..
