#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 1D call option..."
./BSSolverWithStretching solveND cart 1 9 BStest1d.cart.bound BStest1d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest1d.logstretch > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Log Stretching!"
echo "Analytical Value of identical option:"
cat analyticBS.gnuplot | grep "1 "
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 1D call option..."
./BSSolverWithStretching solveND cart 1 9 BStest1d.cart.bound BStest1d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest1d.sinhstretch > tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Sinh Stretching!"
echo "Analytical Value of identical option:"
cat analyticBS.gnuplot | grep "1 "
rm tmp
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 1D put option..."
./BSSolverWithStretching solveND cart 1 9 BStest1d.cart.bound BStest1d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest1d.logstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Log Stretching!"
echo "Analytical Value of identical option:"
cat analyticBS.gnuplot | grep "1 "
rm tmp
echo ""
echo "====================================================================================="
echo "====================================================================================="
echo ""
echo "Executing BSSolverWithStretching: stretched sparse grid with 9 levels, cart, 1D put option..."
./BSSolverWithStretching solveND cart 1 9 BStest1d.cart.bound BStest1d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 2000 0.00001 analytic BStest1d.sinhstretch> tmp
cat tmp | grep "Optionprice at testpoint"
echo "Finished executing BSSolverWithStretching with Sinh Stretching!"
echo "Analytical Value of identical option:"
cat analyticBS.gnuplot | grep "1 "
rm tmp
echo ""
echo "====================================================================================="
echo ""
