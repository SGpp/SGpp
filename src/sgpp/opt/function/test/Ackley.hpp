#ifndef SGPP_OPT_FUNCTION_TEST_ACKLEY_HPP
#define SGPP_OPT_FUNCTION_TEST_ACKLEY_HPP

#include "opt/function/TestFunction.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Ackley : public TestFunction
{
public:
    Ackley(size_t d) : TestFunction(d)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double result = 20.0 + M_E;
        
        double arg1 = 0.0;
        double arg2 = 0.0;
        
        for (size_t t = 0; t < d; t++)
        {
            double xt = 10.0 * x[t] - 1.0;
            arg1 += xt*xt;
            arg2 += cos(2.0 * M_PI * xt);
        }
        
        result += -20.0 * std::exp(-0.2 * std::sqrt(arg1 / static_cast<double>(d)));
        result += -std::exp(arg2 / static_cast<double>(d));
        
        return result;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(d, 0.1);
        return 0.0;
    }
};

}
}
}
}

#endif
