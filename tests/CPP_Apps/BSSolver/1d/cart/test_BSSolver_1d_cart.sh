#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 9 levels, cart, 1D call option..."
./BSSolver solveND cart 1 9 BStest1d.cart.bound BStest1d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint (1.0)"
echo "Finished executing BSSolver!"
echo "Analytical Value of identical option:"
cat analyticBS.gnuplot | grep "1.0"
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolver: regular sparse grid with 9 levels, cart, 1D put option..."
./BSSolver solveND cart 1 9 BStest1d.cart.bound BStest1d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 > tmp
cat tmp | grep "Optionprice at testpoint (1.0)"
echo "Finished executing BSSolver!"
echo "Analytical Value of identical option:"
cat analyticBS.gnuplot | grep "1.0"
rm tmp
echo ""
echo "====================================================================================="
echo ""
