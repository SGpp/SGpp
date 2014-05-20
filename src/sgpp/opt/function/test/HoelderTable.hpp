#ifndef SGPP_OPT_FUNCTION_TEST_HOELDERTABLE_HPP
#define SGPP_OPT_FUNCTION_TEST_HOELDERTABLE_HPP

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

class HoelderTable : public Test
{
public:
    HoelderTable() : Test(2)
    {
    }
    
    void generateDisplacement(unsigned int seed, double standard_deviation)
    {
        Test::generateDisplacement(seed, standard_deviation);
        
        while ((displacement[0] > 0.005) || (displacement[0] < -0.005) ||
               (displacement[1] > 0.01) || (displacement[1] < -0.01))
        {
            seed = std::random_device()();
            Test::generateDisplacement(seed, standard_deviation);
        }
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 20.0 * x[0] - 10.0;
        double x2 = 20.0 * x[1] - 10.0;
        
        return -std::abs(std::sin(x1) * std::cos(x2) *
                         std::exp(std::abs(1.0 - std::sqrt(x1*x1 + x2*x2) / M_PI)));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {0.902751, 0.9832295};
        //return -19.2085;
        return evalUndisplaced(x);
    }
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new HoelderTable(*this));
    }
};

}
}
}
}

#endif
