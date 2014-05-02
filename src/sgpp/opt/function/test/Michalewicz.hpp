#ifndef SGPP_OPT_FUNCTION_TEST_MICHALEWICZ_HPP
#define SGPP_OPT_FUNCTION_TEST_MICHALEWICZ_HPP

#include "opt/function/TestFunction.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Michalewicz : public TestFunction
{
public:
    Michalewicz() : TestFunction(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 5.0 * x[0];
        double x2 = 5.0 * x[1];
        
        return -sin(x1) * std::pow(std::sin(x1*x1 / M_PI), 20.0) -
                sin(x2) * std::pow(std::sin(2.0*x2*x2 / M_PI), 20.0);
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.440638, 0.314098};
        //return -1.8013;
        return evalUndisplaced(x);
    }
};

}
}
}
}

#endif
