/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_LINEAR_NOBOUNDARY_LINEARBASIS_HPP
#define SGPP_OPT_BASIS_LINEAR_NOBOUNDARY_LINEARBASIS_HPP

#include <cmath>
#include <algorithm>

namespace sg
{
namespace opt
{

/**
 * Linear basis on Noboundary grids.
 */
template <class LT, class IT>
class LinearBasis
{
public:
    /**
     * @param l     level of basis function
     * @param i     index of basis function
     * @param x     evaluation point
     * @return      value of linear basis function
     */
    inline double eval(LT l, IT i, double x)
    {
        return std::max(1.0 - std::abs(static_cast<double>(1 << l) * x -
                                       static_cast<double>(i)), 0.0);
    }
};

/// typedef for standard level/index types
typedef LinearBasis<unsigned int, unsigned int> SLinearBase;

}
}

#endif
