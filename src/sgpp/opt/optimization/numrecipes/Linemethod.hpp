#ifndef SGPP_OPT_OPTIMIZATION_NUMRECIPES_LINEMETHOD_HPP
#define SGPP_OPT_OPTIMIZATION_NUMRECIPES_LINEMETHOD_HPP

#include "opt/optimization/numrecipes/Brent.hpp"
#include "opt/optimization/numrecipes/F1dim.hpp"

#include <vector>
#include <cstddef>
#include <cmath>
#include <numeric>

namespace sg
{
namespace opt
{
namespace optimization
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
