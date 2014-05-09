#ifndef SGPP_OPT_FUNCTION_TEST_HARTMAN6_HPP
#define SGPP_OPT_FUNCTION_TEST_HARTMAN6_HPP

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

class Hartman6 : public Test
{
public:
    Hartman6() : Test(6)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        return -1.0 * exp(-10.0*(x[0]-0.1312)*(x[0]-0.1312) - 3.0*(x[1]-0.1696)*(x[1]-0.1696) -
                          17.0*(x[2]-0.5569)*(x[2]-0.5569) - 3.5*(x[3]-0.0124)*(x[3]-0.0124) -
                          1.7*(x[4]-0.8283)*(x[4]-0.8283) - 8.0*(x[5]-0.5886)*(x[5]-0.5886)) -
                1.2 * exp(-0.05*(x[0]-0.2329)*(x[0]-0.2329) - 10.0*(x[1]-0.4135)*(x[1]-0.4135) -
                          17.0*(x[2]-0.8307)*(x[2]-0.8307) - 0.1*(x[3]-0.3736)*(x[3]-0.3736) -
                          8.0*(x[4]-0.1004)*(x[4]-0.1004) - 14.0*(x[5]-0.9991)*(x[5]-0.9991)) -
                3.0 * exp(-3.0*(x[0]-0.2348)*(x[0]-0.2348) - 3.5*(x[1]-0.1451)*(x[1]-0.1451) -
                          1.7*(x[2]-0.3522)*(x[2]-0.3522) - 10.0*(x[3]-0.2883)*(x[3]-0.2883) -
                          17.0*(x[4]-0.3047)*(x[4]-0.3047) - 8.0*(x[5]-0.6650)*(x[5]-0.6650)) -
                3.2 * exp(-17.0*(x[0]-0.4047)*(x[0]-0.4047) - 8.0*(x[1]-0.8828)*(x[1]-0.8828) -
                          0.05*(x[2]-0.8732)*(x[2]-0.8732) - 10.0*(x[3]-0.5743)*(x[3]-0.5743) -
                          0.1*(x[4]-0.1091)*(x[4]-0.1091) - 14.0*(x[5]-0.0381)*(x[5]-0.0381));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573};
        //return -3.322368;
        return evalUndisplaced(x);
    }
};

}
}
}
}

#endif
