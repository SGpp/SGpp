#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 2D call basket option..."
./BSSolverWithStretching solveND cart 2 9 BStest2d.cart.bound BStest2d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest2d.logstretch > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Log Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_call.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 2D call basket option..."
./BSSolverWithStretching solveND cart 2 9 BStest2d.cart.bound BStest2d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest2d.sinhstretch > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Sinh Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_call.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 2D put basket option..."
./BSSolverWithStretching solveND cart 2 9 BStest2d.cart.bound BStest2d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest2d.logstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Log Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_put.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 2D put basket option..."
./BSSolverWithStretching solveND cart 2 9 BStest2d.cart.bound BStest2d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest2d.sinhstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Sinh Stretching!"
echo "Monte Carlo Value of identical option:"
cat mc_value_put.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
