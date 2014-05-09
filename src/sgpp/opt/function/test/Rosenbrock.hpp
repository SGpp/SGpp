#ifndef SGPP_OPT_FUNCTION_TEST_ROSENBROCK_HPP
#define SGPP_OPT_FUNCTION_TEST_ROSENBROCK_HPP

#include "opt/function/test/Test.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Rosenbrock : public Test
{
public:
    Rosenbrock(size_t d) : Test(d)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double result = 0.0;
        
        double xt = 15.0 * x[0] - 5.0;
        double xtm1;
        
        for (size_t t = 1; t < d; t++)
        {
            xtm1 = xt;
            xt = 15.0 * x[t] - 5.0;
            
            double tmp1 = xt - xtm1 * xtm1;
            double tmp2 = 1.0 - xtm1;
            result += 100.0 * tmp1*tmp1 + tmp2*tmp2;
        }
        
        return result;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(d, 0.4);
        return 0.0;
    }
};

}
}
}
}

#endif
