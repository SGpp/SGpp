/**
 * Note by Julian Valentin:
 * 
 * This code is (C) 2007 by Numerical Recipes Software.
 * It is not allowed to be distributed publicly in source code form.
 * Furthermore, it's not allowed to make shared libraries (e.g. DLLs) containing this code
 * or to make commercial binaries containing this code.
 * See http://www.nr.com/licenses/.
 */

#ifndef SGPP_OPT_OPTIMIZER_NUMRECIPES_LINEMETHOD_HPP
#define SGPP_OPT_OPTIMIZER_NUMRECIPES_LINEMETHOD_HPP

#include "opt/optimizer/numrecipes/Brent.hpp"
#include "opt/optimizer/numrecipes/F1dim.hpp"

#include <vector>
#include <cstddef>
#include <cmath>
#include <numeric>

namespace sg
{
namespace opt
{
namespace optimizer
{
namespace numrecipes
{

template <class T>
struct Linemethod
{
    std::vector<double> p, xi;
    T &func;
    size_t n;
    
    Linemethod(T &funcc) : func(funcc)
    {
    }
    
    double linmin()
    {
        double ax, xx, xmin;
        
        n = p.size();
        
        double xi_norm = sqrt(std::inner_product(xi.begin(), xi.end(), xi.begin(), 0.0));
        double d = 0.0;
        
        for (size_t j = 0; j < n; j++)
        {
            double d1 = std::abs(std::min(p[j], 1.0 - p[j]));
            
            if ((j == 0) || (d1 < d))
            {
                d = d1;
            }
        }
        
        F1dim<T> f1dim(p, xi, func);
        ax = 0.0;
        //xx = 1.0;
        //xx = 1.0e-3;
        xx = d / xi_norm;
        
        Brent brent;
        brent.bracket(ax, xx, f1dim);
        xmin = brent.minimize(f1dim);
        
        for (size_t j = 0; j < n; j++)
        {
            xi[j] *= xmin;
            p[j] += xi[j];
        }
        
        return brent.fmin;
    }
};

}
}
}
}

#endif
