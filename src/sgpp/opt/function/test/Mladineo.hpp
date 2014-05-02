#ifndef SGPP_OPT_FUNCTION_TEST_MLADINEO_HPP
#define SGPP_OPT_FUNCTION_TEST_MLADINEO_HPP

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

class Mladineo : public TestFunction
{
public:
    Mladineo() : TestFunction(2)
    {
    }
    
    void generateDisplacement(unsigned int seed, double standard_deviation)
    {
        TestFunction::generateDisplacement(seed, standard_deviation);
        
        while ((displacement[0] > 0) || (displacement[0] < -0.01) ||
               (displacement[1] > 0) || (displacement[1] < -0.01))
        {
            seed = std::random_device()();
            TestFunction::generateDisplacement(seed, standard_deviation);
        }
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 0.99 * x[0] + 0.01;
        double x2 = 0.99 * x[1] + 0.01;
        
        return 1.0 + (x1*x1 + x2*x2) / 2.0 - cos(10.0*log(2.0*x1))*cos(10.0*log(3.0*x2));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.001542, 0.004449};
        //return 0.000170;
        return evalUndisplaced(x);
    }
};

}
}
}
}

#endif
