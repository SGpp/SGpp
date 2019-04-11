
## Monomials
# time python2.7 comparingBSplines.py --model=monomial --scalarModelParameter=3 --dim=1 --gridType=nak --degree=135 --refineType=regularByPoints --maxPoints=1000
# time python2.7 comparingBSplines.py --model=monomial --scalarModelParameter=3 --dim=1 --gridType=nak --degree=135 --refineType=surplus --initialLevel=1 --maxPoints=1000 --numSteps=7

## plot all monomials
# python2.7 plotBsplineComparison.py --scalarModelParameter=0
# python2.7 plotBsplineComparison.py --scalarModelParameter=1
# python2.7 plotBsplineComparison.py --scalarModelParameter=2
# python2.7 plotBsplineComparison.py --scalarModelParameter=3
# python2.7 plotBsplineComparison.py --scalarModelParameter=4
# python2.7 plotBsplineComparison.py --scalarModelParameter=5

## Exp
# time python2.7 comparingBSplines.py --model=sumExp --dim=10 --gridType=naknobound --degree=35 --refineType=regularByPoints --initialLevel=2 --maxPoints=10000
# time python2.7 comparingBSplines.py --model=sumExp --dim=10 --gridType=naknobound --degree=35 --refineType=surplus --initialLevel=2 --maxPoints=10000 --numSteps=7

## Sin
# time python2.7 comparingBSplines.py --model=prodSin --dim=4 --gridType=naknobound --degree=135 --refineType=regularByPoints --initialLevel=2 --maxPoints=4000
# time python2.7 comparingBSplines.py --model=prodSin --dim=4 --gridType=naknobound --degree=135 --refineType=surplus --initialLevel=2 --maxPoints=10000 --numSteps=7

## Cos
# time python2.7 comparingBSplines.py --model=cosSum --dim=10 --gridType=naknobound --degree=135 --refineType=regularByPoints --initialLevel=2 --maxPoints=10000
# time python2.7 comparingBSplines.py --model=cosSum --dim=10 --gridType=naknobound --degree=135 --refineType=surplus --initialLevel=2 --maxPoints=10000 --numSteps=7

## wing
#time python2.7 comparingBSplines.py --model=wing --gridType=nakbsplineextended --degree=3 --refineType=surplus --initialLevel=2 --maxPoints=1000 --numSteps=5 --numRefine=20 --saveData=0

time python2.7 comparingBSplines.py --model=dette --gridType=nak --degree=135 --refineType=surplus --initialLevel=2 --maxPoints=10000 --numSteps=7
time python2.7 comparingBSplines.py --model=dette --gridType=nak --degree=135 --refineType=regularByPoints --maxPoints=10000 

## borehole
# time python2.7 comparingBSplines.py --model=borehole --gridType=nak --degree=135 --refineType=regularByPoints --maxPoints=10000 
# time python2.7 comparingBSplines.py --model=borehole --gridType=nak --degree=135 --refineType=surplus --initialLevel=2 --maxPoints=10000 --numSteps=7






