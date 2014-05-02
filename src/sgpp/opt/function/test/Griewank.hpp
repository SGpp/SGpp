#ifndef SGPP_OPT_FUNCTION_TEST_GRIEWANK_HPP
#define SGPP_OPT_FUNCTION_TEST_GRIEWANK_HPP

#include "opt/function/TestFunction.hpp"

#include <cmath>

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Griewank : public TestFunction
{
public:
    Griewank(size_t d) : TestFunction(d)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double result = 1.0;
        double tmp = 1.0;
        
        for (size_t t = 0; t < d; t++)
        {
            double xt = 1200.0 * x[t] - 600.0;
            result += xt*xt / 4000.0;
            tmp *= cos(xt / sqrt(static_cast<double>(t+1)));
        }
        
        result -= tmp;
        return result;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(d, 0.5);
        return 0.0;
    }
};

}
}
}
}

#endif
