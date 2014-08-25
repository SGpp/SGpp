/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_LINEAR_BOUNDARY_LINEARBOUNDARYBASIS_HPP
#define SGPP_OPT_BASIS_LINEAR_BOUNDARY_LINEARBOUNDARYBASIS_HPP

#include <cmath>

namespace sg
{
namespace opt
{

/**
 * Linear basis on Boundary grids.
 */
template <class LT, class IT>
class LinearBoundaryBasis
{
public:
    /**
     * @param l     level of basis function
     * @param i     index of basis function
     * @param x     evaluation point
     * @return      value of boundary linear basis function
     */
    inline double eval(LT l, IT i, double x)
    {
        if (l == 0)
        {
            // first level
            if (i == 0)
            {
                return 1.0 - x;
            } else
            {
                return x;
            }
        } else
        {
            return std::max(1.0 - std::abs(static_cast<double>(1 << l) * x -
                                           static_cast<double>(i)), 0.0);
        }
    }
};

/// typedef for standard level/index types
typedef LinearBoundaryBasis<unsigned int, unsigned int> SLinearBoundaryBase;

}
}

#endif
