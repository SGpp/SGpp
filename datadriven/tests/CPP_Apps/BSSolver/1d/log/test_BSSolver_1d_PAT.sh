#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 9 levels, PAT, 1D call option..."
./BSSolver solveND PAT 1 9 BStest1d.log.bound BStest1d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolver!"
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 9 levels, PAT, 1D put option..."
./BSSolver solveND PAT 1 9 BStest1d.log.bound BStest1d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolver!"
rm tmp
echo ""
echo "====================================================================================="
echo ""
