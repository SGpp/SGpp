#!/bin/sh
echo ""
echo "====================================================================================="
echo ""
echo "Executing BSHWSolver: regular sparse grid with 5 levels, cart, call payoff (Black-Scholes / Hull-White)"
./BSHWSolver CombineBSHW 2 5 0.01 0.1 BSHW2d.stoch.data BSHW2d.bound.data std_euro_call 1 0.015625 400 0.000001 ImEul 1.0 0.3 cart > tmp
cat tmp | grep "Optionprice"
echo "Finished executing BSHWSolver!"
echo "Reference value: "
cat ref_value.txt
rm tmp
echo ""
echo "====================================================================================="
echo ""
