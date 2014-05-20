#ifndef SGPP_OPT_FUNCTION_TEST_EGGHOLDER_HPP
#define SGPP_OPT_FUNCTION_TEST_EGGHOLDER_HPP

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

class Eggholder : public Test
{
public:
    Eggholder() : Test(2)
    {
    }
    
    void generateDisplacement(unsigned int seed, double standard_deviation)
    {
        Test::generateDisplacement(seed, standard_deviation);
        displacement[0] = 0.0;
    }
    
    double evalUndisplaced(const std::vector<double> &x)
    {
        double x1 = 1024.0 * x[0] - 512.0;
        double x2 = 1024.0 * x[1] - 512.0;
        
        return -(x2 + 47.0) * std::sin(std::sqrt(std::abs(x1/2.0 + x2 + 47.0))) -
                x1 * std::sin(std::sqrt(std::abs(x1 - (x2 + 47.0))));
    }
    
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = {1.0, 0.8947577};
        //return -959.6407;
        return evalUndisplaced(x);
    }
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new Eggholder(*this));
    }
};

}
}
}
}

#endif
