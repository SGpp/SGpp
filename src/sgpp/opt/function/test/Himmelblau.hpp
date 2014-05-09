#ifndef SGPP_OPT_FUNCTION_TEST_HIMMELBLAU_HPP
#define SGPP_OPT_FUNCTION_TEST_HIMMELBLAU_HPP

#include "opt/function/test/Test.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Himmelblau : public Test
{
public:
    Himmelblau() : Test(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 10.0 * x[0] - 5.0;
        double x2 = 10.0 * x[1] - 5.0;
        
        return (x1*x1 + x2 - 11.0)*(x1*x1 + x2 - 11.0) + (x1 + x2*x2 - 7.0)*(x1 + x2*x2 - 7.0);
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.8, 0.7};
        return 0.0;
    }
};

}
}
}
}

#endif
