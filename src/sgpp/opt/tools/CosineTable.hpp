/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_TOOLS_COSINETABLE_HPP
#define SGPP_OPT_TOOLS_COSINETABLE_HPP

#include <vector>

namespace sg
{
namespace opt
{
namespace tools
{

/**
 * Lookup table for cosine values.
 * When working with Clenshaw-Curtis grids, cosine values are needed frequently.
 * This class divides \f$[0, \pi/2]\f$ uniformly into \f$n\f$ intervals,
 * evaluates on construction the cosine function at the interval endpoints
 * and interpolates linearly when doing lookup operations.
 */
class CosineTable
{
public:
    /// default number of intervals
    static const size_t DEFAULT_SIZE = 1 << 20;
    
    /**
     * Constructor creating the lookup table.
     * 
     * @param size      number of sub-intervals of \f$[0, \pi/2]\f$
     */
    CosineTable(size_t size = DEFAULT_SIZE) :
        table(std::vector<double>(size+1, 0.0)),
        h(M_PI_2 / static_cast<double>(size)),
        hinv(1.0 / h)
    {
        double x = 0.0;
        
        for (size_t i = 0; i <= size; i++)
        {
            table[i] = cos(x);
            // rounding errors can be neglected
            x += h;
        }
    }
    
    /**
     * @param x     cosine argument, must be in \f$[0, \pi]\f$
     */
    inline double lookUp(double x) const
    {
        double factor = 1.0;
        
        if (x > M_PI_2)
        {
            // mirror the situation at x = pi/2
            x = M_PI - x;
            factor = -1.0;
        } else if (x == M_PI_2)
        {
            return 0.0;
        }
        
        const double xl = x * hinv;
        // x should be in [table[i], table[i+1])
        const size_t i = static_cast<size_t>(xl);
        // convex combination coefficient
        const double t = xl - static_cast<double>(i);
        
        // linear interpolation
        return factor * ((1.0-t) * table[i] + t * table[i+1]);
    }
    
protected:
    /// lookup table
    std::vector<double> table;
    /// interval length
    double h;
    /// \f$1/h\f$
    double hinv;
};

}
}
}

#endif
