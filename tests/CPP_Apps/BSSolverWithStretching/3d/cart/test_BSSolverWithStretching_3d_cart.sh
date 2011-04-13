#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 8 levels, cart, 3D call basket option..."
./BSSolverWithStretching solveND cart 3 8 BStest3d.cart.bound BStest3d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest3d.logstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Log Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_call.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 8 levels, cart, 3D call basket option..."
./BSSolverWithStretching solveND cart 3 8 BStest3d.cart.bound BStest3d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest3d.sinhstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Sinh Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_call.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 8 levels, cart, 3D put basket option..."
./BSSolverWithStretching solveND cart 3 8 BStest3d.cart.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest3d.logstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Log Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_put.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 8 levels, cart, 3D put basket option..."
./BSSolverWithStretching solveND cart 3 8 BStest3d.cart.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest3d.sinhstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Sinh Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_put.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
