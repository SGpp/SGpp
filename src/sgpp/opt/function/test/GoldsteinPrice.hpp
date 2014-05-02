#ifndef SGPP_OPT_FUNCTION_TEST_GOLDSTEINPRICE_HPP
#define SGPP_OPT_FUNCTION_TEST_GOLDSTEINPRICE_HPP

#include "opt/function/TestFunction.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class GoldsteinPrice : public TestFunction
{
public:
    GoldsteinPrice() : TestFunction(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 4.0 * x[0] - 2.0;
        double x2 = 4.0 * x[1] - 2.0;
        
        return (1.0 + (x1+x2+1.0)*(x1+x2+1.0) *
                      (19.0 - 14.0*x1 + 3.0*x1*x1 - 14.0*x2 + 6.0*x1*x2 + 3.0*x2*x2)) *
               (30.0 + (2.0*x1-3.0*x2)*(2.0*x1-3.0*x2) *
                       (18.0 - 32.0*x1 + 12.0*x1*x1 + 48.0*x2 - 36.0*x1*x2 + 27.0*x2*x2));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.5, 0.25};
        return 3.0;
    }
};

}
}
}
}

#endif
