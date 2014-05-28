/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Florian Zipperle (florian.zipperle@tum.de)

#include <algorithm>    // std::swap
#include <math.h>
#include "OperationL2ScalarProduct.hpp"
using namespace std;
/*

 def phi_1d(l, i, x):
 """Evalute the point in 1d for periodic basis set"""
 if l == 0:
 if x >= 0 and x < 0.5 :
 return 1.0 - 2 * x
 elif x >= 0.5 and x <= 1:
 return 2 * (x - 0.5)
 else:
 return 0.0
 else:
 return max([1 - abs(2 ** l * x - i), 0])

 */

/*

 def get_ancestor(l1, l2, i2):
 """Calculate the ancestor grid function of (l2,i2) on the level l1"""
 result = int(i2 / (2 ** (l2 - l1)))
 result += 1 - result % 2
 return result

 */

namespace sg {
	namespace datadriven {
		//"""Calculate the ancestor grid function of (l2,i2) on the level l1"""
int get_ancestor(int l1, int l2, int i2)
{
    int result = int(i2 / (pow(2.0, (l2 - l1))));
    result += 1 - result % 2;
    return result;
}

double L2_1d(int l1, int i1, int l2, int i2)
{

    if (l1 > l2)
    {
        swap(l1, l2);
        swap(i1, i2);
    }

    int anc_idx = get_ancestor(max(l1, 1), l2, i2);

    //case 0: same level, different functions
    if (l1 == l2 && i1 != i2)
    {
        return 0.0;
    }
    else if (l1 == 0 && l2 == 0)
    {
        return 1.0 / 3;
    }
    else if (l1 == 0 && l2 == 1)
    {
        return 1.0 / 6;
    }
    else if (l1 == 0 && (i2 == 1 || i2 == (pow(2.0,l2)) - 1))
    { //case 7: border part of folded out function and boundary hat function
        return pow(2.0,-l2 - 1);
    }
    else if (l1 == l2 && i1 == i2)
    { //case 4: twice the same function
        return (pow(2.0,1 - l1)) / 3.0;
    }
    else if (i1 == anc_idx || l1 == 0)
    { //case 5: two hat functions OR a folded out and a hat functions
        bool isLvlZero = false;
        if (l1 == 0)
        {
            isLvlZero = true;
            l1 = 1;
            i1 = 1;
        }

        double result = pow(2.0,l1 - 2 * l2);

        double step = pow(2.0,-(l2 - l1 - 1));

        int i1_int = int(i1 / 2) + 1;
        int i2_int = int(i2 / 2) + 1;

        int i = i2_int - (i1_int - 1) * int(pow(2.0,l2 - l1));
        int mid = int(pow(2.0,l2 - l1 - 1));

        int mult = 1;
        if (i <= mid)
        {
            mult = i - 1;
        }
        else
        {
            mult = mid * 2 - i;
        }
        if (isLvlZero)
        {
            mult = (mid - 1) - mult;
        }
        double A = pow(2.0,-l2);

        result += A * mult * step;
        return result;
    }
    else
    {
        return 0.0;
    }
}
	}
}
/*
 def L2_1d(l1, i1, l2, i2):
 """Evaluate the L2 product in 1d for periodic basis set"""

 if l1 > l2:
 return L2_1d(l2, i2, l1, i1)

 anc_idx = get_ancestor(np.max([l1, 1]), l2, i2)

 # case 0: same level, different functions
 if l1 == l2 and i1 != i2:
 return 0.0

 elif l1 == 0 and l2 == 0:
 return 1.0 / 3

 elif l1 == 0 and l2 == 1:
 return 1.0 / 6

 # case 7: border part of folded out function and boundary hat function
 # elif l1 == 0 and (i2 == 1 or i2 ==2**l2-1 ):
 #    return 2**(-l2-1)

 # case 4: twice the same function
 elif l1 == l2 and i1 == i2:
 return 2 ** (1 - l1) / 3.0

 # case 5: two hat functions OR a folded out and a hat functions
 elif i1 == anc_idx or l1 == 0:
 isLvlZero = False
 if l1 == 0:
 isLvlZero = True
 l1 = 1
 i1 = 1

 result = 2 ** (l1 - 2 * l2)

 step = 2 ** -(l2 - l1 - 1)

 i1_int = int(i1 / 2) + 1
 i2_int = int(i2 / 2) + 1

 i = i2_int - (i1_int - 1) * 2 ** (l2 - l1)
 mid = 2 ** (l2 - l1 - 1)


 if(i <= mid):
 mult = i - 1
 else:
 mult = mid * 2 - i

 if isLvlZero:
 mult = (mid - 1) - mult

 A = 2 ** -l2

 result += A * mult * step
 return result

 else:
 return 0.0
 */
