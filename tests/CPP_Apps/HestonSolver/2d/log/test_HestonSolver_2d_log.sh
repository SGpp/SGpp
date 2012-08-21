#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing HestonSolver: regular sparse grid with 4 levels, log, 1D call option..."
./HestonSolver solveND log 2 4 Hestontest2d.log.bound Hestontest2d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing HestonSolver!"
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing HestonSolver: regular sparse grid with 4 levels, log, 1D put option..."
./HestonSolver solveND log 2 4 Hestontest2d.log.bound Hestontest2d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing HestonSolver!"
rm tmp
echo ""
echo "====================================================================================="
echo ""
