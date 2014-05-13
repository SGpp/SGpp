#ifndef SGPP_OPT_FUNCTION_TEST_BEALE_HPP
#define SGPP_OPT_FUNCTION_TEST_BEALE_HPP

#include "opt/function/test/Test.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class Beale : public Test
{
public:
    Beale() : Test(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 10.0 * x[0] - 5.0;
        double x2 = 10.0 * x[1] - 5.0;
        double tmp1 = 1.5 - x1 * (1.0 - x2);
        double tmp2 = 2.25 - x1 * (1.0 - x2*x2);
        double tmp3 = 2.625 - x1 * (1.0 - x2*x2*x2);
        
        return tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.8, 0.55};
        return 0.0;
    }
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new Beale(*this));
    }
};

}
}
}
}

#endif
