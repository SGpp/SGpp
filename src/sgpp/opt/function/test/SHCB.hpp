#ifndef SGPP_OPT_FUNCTION_TEST_SHCB_HPP
#define SGPP_OPT_FUNCTION_TEST_SHCB_HPP

#include "opt/function/test/Test.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

class SHCB : public Test
{
public:
    SHCB() : Test(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 10.0 * x[0] - 5.0;
        double x2 = 10.0 * x[1] - 5.0;
        
        return x1*x1 * (4.0 - 2.1*x1*x1 + x1*x1*x1*x1/3.0) + x1*x2 + 4.0*x2*x2 * (x2*x2 - 1.0);
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.508984, 0.428734};
        //return -1.031628;
        return evalUndisplaced(x);
    }
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new SHCB(*this));
    }
};

}
}
}
}

#endif
