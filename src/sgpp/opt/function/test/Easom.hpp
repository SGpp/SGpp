#ifndef SGPP_OPT_FUNCTION_TEST_EASOM_HPP
#define SGPP_OPT_FUNCTION_TEST_EASOM_HPP

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

class Easom : public Test
{
public:
    Easom() : Test(2)
    {
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 200.0 * x[0] - 100.0;
        double x2 = 200.0 * x[1] - 100.0;
        
        return -std::cos(x1) * std::cos(x2) *
                std::exp(-((x1-M_PI)*(x1-M_PI) + (x2-M_PI)*(x2-M_PI)));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.51570796326794896619231, 0.51570796326794896619231};
        return -1.0;
    }
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new Easom(*this));
    }
};

}
}
}
}

#endif
