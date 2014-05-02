/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef BSPLINE_BOUNDARY_BASE_HPP
#define BSPLINE_BOUNDARY_BASE_HPP

#include "base/basis/bspline/noboundary/BsplineBasis.hpp"

#include <cmath>

namespace sg
{
namespace base
{

template<class LT, class IT>
class BsplineBoundaryBasis
{
protected:
    BsplineBasis<LT, IT> bspline_basis;
    
public:
    BsplineBoundaryBasis() : bspline_basis(BsplineBasis<LT, IT>()) {}
    BsplineBoundaryBasis(size_t degree) : bspline_basis(BsplineBasis<LT, IT>(degree)) {}
    
    inline double eval(LT level, IT index, double x) const
    {
        if (level == 0)
        {
            return bspline_basis.uniformBSpline(
                    x + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        } else
        {
            double hinv = static_cast<double>(1 << level);
            
            return bspline_basis.uniformBSpline(
                    x * hinv + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        }
    }
    
    inline double evalDx(LT level, IT index, double x) const
    {
        if (level == 0)
        {
            return bspline_basis.uniformBSplineDx(
                    x + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        } else
        {
            double hinv = static_cast<double>(1 << level);
            
            return hinv * bspline_basis.uniformBSplineDx(
                    x * hinv + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        }
    }
    
    inline double evalDxDx(LT level, IT index, double x) const
    {
        if (level == 0)
        {
            return bspline_basis.uniformBSplineDxDx(
                    x + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        } else
        {
            double hinv = static_cast<double>(1 << level);
            
            return hinv*hinv * bspline_basis.uniformBSplineDxDx(
                    x * hinv + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        }
    }
    
    inline size_t getDegree() const { return bspline_basis.degree; }
};

typedef BsplineBoundaryBasis<unsigned int, unsigned int> SBsplineBoundaryBase;

}
}

#endif /* BSPLINE_BOUNDARY_BASE_HPP */
