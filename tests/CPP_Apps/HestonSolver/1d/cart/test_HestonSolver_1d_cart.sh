#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing HestonSolver: regular sparse grid with 9 levels, cart, 1D call option..."
./HestonSolver solveND cart 1 9 Hestontest1d.cart.bound Hestontest1d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing HestonSolver!"
echo "Analytical Value of identical option:"
cat analyticHeston.gnuplot | grep "1 "
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing HestonSolver: regular sparse grid with 9 levels, cart, 1D put option..."
./HestonSolver solveND cart 1 9 Hestontest1d.cart.bound Hestontest1d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing HestonSolver!"
echo "Analytical Value of identical option:"
cat analyticHeston.gnuplot | grep "1 "
rm tmp
echo ""
echo "====================================================================================="
echo ""
