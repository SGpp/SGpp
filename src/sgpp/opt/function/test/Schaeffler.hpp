#ifndef SGPP_OPT_FUNCTION_TEST_SCHAEFFLER_HPP
#define SGPP_OPT_FUNCTION_TEST_SCHAEFFLER_HPP

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

class Schaeffler : public Test
{
public:
    Schaeffler() : Test(50)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 4.0 * x[0] - 1.05;
        double result = 1.0 + 6.0*x1*x1 - cos(12.0*x1);
        
        double xi = x1;
        double xim1;
        
        for (size_t i = 1; i < 50; i++)
        {
            xim1 = xi;
            xi = 4.0 * x[i] - 1.05;
            
            result += 590.0 * (xi-xim1)*(xi-xim1);
        }
        
        return result;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(50, 0.2625);
        return 0.0;
    }
};

}
}
}
}

#endif
