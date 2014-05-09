#ifndef SGPP_OPT_OPTIMIZER_ARMIJORULE_HPP
#define SGPP_OPT_OPTIMIZER_ARMIJORULE_HPP

#include "base/datatypes/DataVector.hpp"
#include "opt/function/Objective.hpp"

#include <cstddef>

namespace sg
{
namespace opt
{
namespace optimizer
{

inline bool armijoRule(function::Objective &f, double beta, double gamma,
                       const std::vector<double> &x, double fx, base::DataVector &grad_fx,
                       const std::vector<double> &s, std::vector<double> &y)
{
    size_t d = x.size();
    double sigma = 1.0;
    double ip = 0.0;
    double fy = fx;
    
    for (size_t t = 0; t < d; t++)
    {
        ip += grad_fx[t] * s[t];
    }
    
    for (size_t k = 0; k < 100; k++)
    {
        bool y_in_domain = true;
        
        for (size_t t = 0; t < d; t++)
        {
            y[t] = x[t] + sigma * s[t];
            
            if ((y[t] < 0.0) || (y[t] > 1.0))
            {
                y_in_domain = false;
                break;
            }
        }
        
        if (y_in_domain)
        {
            fy = f.eval(y);
            
            //std::cout << "sigma: " << sigma << ", fy: " << fy << "\n";
            if (fy - fx <= sigma * gamma * ip)
            {
                /*std::cout << "ip: " << ip << "\n";
                std::cout << "sigma: " << sigma << "\n";
                std::cout << "fy: " << fy << "\n";*/
                return true;
            }
        }
        
        sigma *= beta;
    }
    
    return false;
}

}
}
}

#endif
