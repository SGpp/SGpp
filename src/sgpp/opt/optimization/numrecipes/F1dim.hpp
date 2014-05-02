#ifndef SGPP_OPT_OPTIMIZATION_NUMRECIPES_F1DIM_HPP
#define SGPP_OPT_OPTIMIZATION_NUMRECIPES_F1DIM_HPP

#include <vector>
#include <cstddef>
#include <cmath>

namespace sg
{
namespace opt
{
namespace optimization
{
namespace numrecipes
{

template <class T>
struct F1dim
{
    const std::vector<double> &p;
    const std::vector<double> &xi;
    //std::vector<double> xi;
    size_t n;
    T &func;
    std::vector<double> xt;
    
    F1dim(const std::vector<double> &pp, const std::vector<double> &xii, T &funcc) :
            p(pp), xi(xii), n(pp.size()), func(funcc), xt(n)
            //p(pp), n(pp.size()), func(funcc), xt(n)
    {
        /*xi_orig_norm = sqrt(std::inner_product(xii.begin(), xii.end(), xii.begin(), 0.0));
        xi = std::vector<double>(n, 0.0);
        
        for (int j = 0; j < n; j++)
        {
            xi[j] = xii[j] / xi_orig_norm;
        }*/
    }
    
    double operator()(const double x)
    {
        for (size_t j = 0; j < n; j++)
        {
            xt[j] = p[j] + x * xi[j];
        }
        
        return func(xt);
    }
};

}
}
}
}

#endif
