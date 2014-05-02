#ifndef SGPP_OPT_FUNCTION_TEST_BRANIN_HPP
#define SGPP_OPT_FUNCTION_TEST_BRANIN_HPP

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

class Branin : public TestFunction
{
public:
    Branin() : TestFunction(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 15.0 * x[0] - 5.0;
        double x2 = 15.0 * x[1];
        double tmp = x2 - 5.1 * x1 * x1 / (4.0 * M_PI * M_PI) + 5.0 * x1 / M_PI - 6.0;
        
        return tmp * tmp + 10.0 * (1.0 - 1.0 / (8.0 * M_PI)) * std::cos(x1) + 10.0;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.961652, 0.165};
        //return 0.397887;
        return evalUndisplaced(x);
    }
};

}
}
}
}

#endif
