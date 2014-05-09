/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef WAVELET_BASE_HPP
#define WAVELET_BASE_HPP

#include "base/basis/Basis.hpp"

#include <cmath>

namespace sg
{
namespace base
{

template<class LT, class IT>
class WaveletBasis : public Basis<LT, IT>
{
public:
    inline double eval(LT level, IT index, double x)
    {
        double hinv = static_cast<double>(1 << level);
        double t = x * hinv - static_cast<double>(index);
        
        if ((t > 2.0) || (t < -2.0))
        {
            return 0.0;
        }
        
        double t2 = t*t;
        
        return (1.0 - t2) * exp(-t2);
    }
    
    inline double evalDx(LT level, IT index, double x)
    {
        double hinv = static_cast<double>(1 << level);
        double t = x * hinv - static_cast<double>(index);
        
        if ((t > 2.0) || (t < -2.0))
        {
            return 0.0;
        }
        
        double t2 = t*t;
        
        return 2.0 * t * (t2 - 2.0) * exp(-t2) * hinv;
    }
    
    inline double evalDxDx(LT level, IT index, double x)
    {
        double hinv = static_cast<double>(1 << level);
        double t = x * hinv - static_cast<double>(index);
        
        if ((t > 2.0) || (t < -2.0))
        {
            return 0.0;
        }
        
        double t2 = t*t;
        double t4 = t2*t2;
        
        return -2.0 * (2.0*t4 - 7.0*t2 + 2.0) * exp(-t2) * hinv*hinv;
    }
};

typedef WaveletBasis<unsigned int, unsigned int> SWaveletBase;

}
}

#endif /* WAVELET_BASE_HPP */
