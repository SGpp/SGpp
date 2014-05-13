/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef WAVELET_BOUNDARY_BASE_HPP
#define WAVELET_BOUNDARY_BASE_HPP

#include "base/basis/Basis.hpp"

#include <cmath>

namespace sg
{
namespace base
{

template<class LT, class IT>
class WaveletBoundaryBasis : public Basis<LT, IT>
{
public:
    inline double eval(LT level, IT index, double x)
    {
        double hinv = static_cast<double>(1 << level);
        double t = x * hinv - static_cast<double>(index);
        
        if ((t >= 2.0) || (t <= -2.0))
        {
            return 0.0;
        }
        
        double t2 = t*t;
        
        return (1.0 - t2) * exp(-t2);
        //return (1.0 - t2) * exp(-t2 + 1.0/(t2-4.0) + 0.25);
    }
    
    inline double evalDx(LT level, IT index, double x)
    {
        double hinv = static_cast<double>(1 << level);
        double t = x * hinv - static_cast<double>(index);
        
        if ((t >= 2.0) || (t <= -2.0))
        {
            return 0.0;
        }
        
        double t2 = t*t;
        //double s = 1.0/(t2-4.0);
        
        return 2.0 * t * (t2 - 2.0) * exp(-t2) * hinv;
        //return 2.0 * t * ((t2-1.0) * (s*s + 1.0) - 1.0) * exp(-t2 + s + 0.25) * hinv;
    }
    
    inline double evalDxDx(LT level, IT index, double x)
    {
        double hinv = static_cast<double>(1 << level);
        double t = x * hinv - static_cast<double>(index);
        
        if ((t >= 2.0) || (t <= -2.0))
        {
            return 0.0;
        }
        
        double t2 = t*t;
        //double s = 1.0/(t2-4.0);
        
        return -2.0 * (2.0*t2*t2 - 7.0*t2 + 2.0) * exp(-t2) * hinv*hinv;
        /*return -2.0 * (2.0 * (t2-1.0) * t2 * (s*s + 1.0) * (s*s + 1.0)
                       - 4.0 * t2 * (s*s + 1.0)
                       + (t2-1.0) * (4.0*t2*s*s*s - s*s - 1.0)
                       + 1.0)
                * exp(-t2 + s + 0.25) * hinv*hinv;*/
    }
};

typedef WaveletBoundaryBasis<unsigned int, unsigned int> SWaveletBoundaryBase;

}
}

#endif /* WAVELET_BOUNDARY_BASE_HPP */
