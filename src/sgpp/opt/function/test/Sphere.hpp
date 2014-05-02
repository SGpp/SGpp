#ifndef SGPP_OPT_FUNCTION_TEST_SPHERE_HPP
#define SGPP_OPT_FUNCTION_TEST_SPHERE_HPP

#include "opt/function/TestFunction.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Sphere : public TestFunction
{
public:
    Sphere(size_t d) : TestFunction(d)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double result = 0.0;
        
        for (size_t t = 0; t < d; t++)
        {
            double xt = 10.24 * x[t] - 5.12;
            result += xt*xt;
        }
        
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
