/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_LINEAR_MODIFIED_MODIFIEDLINEARBASIS_HPP
#define SGPP_OPT_BASIS_LINEAR_MODIFIED_MODIFIEDLINEARBASIS_HPP

#include <cmath>

namespace sg
{
namespace opt
{

/**
 * Modified linear basis on Noboundary grids.
 */
template <class LT, class IT>
class ModLinearBasis
{
public:
    /**
     * @param l     level of basis function
     * @param i     index of basis function
     * @param x     evaluation point
     * @return      value of modified linear basis function
     */
    inline double eval(LT l, IT i, double x)
    {
        const IT hinv = 1 << l;
        const double h = 1.0 / static_cast<double>(hinv);
        
        if (l == 1)
        {
            // first level
            return 1.0;
        } else if (i == 1)
        {
            // left modified basis function
            return ((x <= 2.0 * h) ? (2.0 - hinv * x) : 0.0);
        } else if (i == hinv - 1)
        {
            // right modified basis function
            return ((x >= 1.0 - 2.0 * h) ? (hinv * x - i + 1.0) : 0.0);
        } else
        {
            // interior basis function
            return std::max(1.0 - std::abs(hinv * x - i), 0.0);
        }
    }
};

/// typedef for standard level/index types
typedef ModLinearBasis<unsigned int, unsigned int> SModLinearBase;

}
}

#endif

