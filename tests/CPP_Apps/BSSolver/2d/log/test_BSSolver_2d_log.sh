#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 9 levels, log, 2D call basket option..."
./BSSolver solveND log 2 9 BStest2d.log.bound BStest2d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolver!"
echo "Monte Carlo Value of identical option:"
cat mc_value_call.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 9 levels, log, 2D put basket option..."
./BSSolver solveND log 2 9 BStest2d.log.bound BStest2d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolver!"
echo "Monte Carlo Value of identical option:"
cat mc_value_put.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
