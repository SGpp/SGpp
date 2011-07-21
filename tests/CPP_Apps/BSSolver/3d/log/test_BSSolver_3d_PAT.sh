#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 8 levels, PAT, 3D call basket option..."
./BSSolver solveND PAT 3 8 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolver!"
echo "Monte Carlo Value of identical option:"
cat mc_value_call.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 8 levels, PAT, 3D put basket option..."
./BSSolver solveND PAT 3 8 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolver!"
echo "Monte Carlo Value of identical option:"
cat mc_value_put.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
