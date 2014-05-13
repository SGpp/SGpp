#ifndef SGPP_OPT_FUNCTION_TEST_RASTRIGIN_HPP
#define SGPP_OPT_FUNCTION_TEST_RASTRIGIN_HPP

#include "opt/function/test/Test.hpp"

#include <cmath>

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Rastrigin : public Test
{
public:
    Rastrigin(size_t d) : Test(d)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double result = 10.0 * static_cast<double>(d);
        
        for (size_t t = 0; t < d; t++)
        {
            double xt = 10.0 * x[t] - 2.0;
            result += xt*xt - 10.0 * cos(2*M_PI*xt);
        }
        
        return result;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(d, 0.2);
        return 0.0;
    }
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new Rastrigin(*this));
    }
};

}
}
}
}

#endif
