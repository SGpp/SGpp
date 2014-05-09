#ifndef SGPP_OPT_FUNCTION_TEST_HARTMAN3_HPP
#define SGPP_OPT_FUNCTION_TEST_HARTMAN3_HPP

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

class Hartman3 : public Test
{
public:
    Hartman3() : Test(3)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        return -1.0 * exp(-3.0*(x[0]-0.3689)*(x[0]-0.3689) - 10.0*(x[1]-0.1170)*(x[1]-0.1170) -
                          30.0*(x[2]-0.2673)*(x[2]-0.2673)) -
                1.2 * exp(-0.1*(x[0]-0.4699)*(x[0]-0.4699) - 10.0*(x[1]-0.4387)*(x[1]-0.4387) -
                          35.0*(x[2]-0.7470)*(x[2]-0.7470)) -
                3.0 * exp(-3.0*(x[0]-0.1091)*(x[0]-0.1091) - 10.0*(x[1]-0.8732)*(x[1]-0.8732) -
                          30.0*(x[2]-0.5547)*(x[2]-0.5547)) -
                3.2 * exp(-0.1*(x[0]-0.0382)*(x[0]-0.0382) - 10.0*(x[1]-0.5743)*(x[1]-0.5743) -
                          35.0*(x[2]-0.8828)*(x[2]-0.8828));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.114614, 0.555649, 0.852547};
        //return -3.862785;
        return evalUndisplaced(x);
    }
};

}
}
}
}

#endif
