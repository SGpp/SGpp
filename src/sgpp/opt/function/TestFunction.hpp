#ifndef SGPP_OPT_FUNCTION_TESTFUNCTION_HPP
#define SGPP_OPT_FUNCTION_TESTFUNCTION_HPP

#include "opt/function/ObjectiveFunction.hpp"

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace function
{

class TestFunction : public ObjectiveFunction
{
public:
    static const double DEFAULT_STANDARD_DEVIATION;
    
    TestFunction(size_t d);
    virtual ~TestFunction() {}
    
    double eval(const std::vector<double> &x);
    virtual double evalUndisplaced(const std::vector<double> &x) = 0;
    
    double getOptimalPoint(std::vector<double> &x);
    virtual double getOptimalPointUndisplaced(std::vector<double> &x) = 0;
    
    void generateDisplacement();
    void generateDisplacement(unsigned int seed);
    virtual void generateDisplacement(unsigned int seed, double standard_deviation);
    
    void displaceVector(std::vector<double> &x) const;
    void reverseDisplaceVector(std::vector<double> &x) const;
    
    unsigned int getSeed() const;
    double getStandardDeviation() const;
    void getDisplacement(std::vector<double> &displacement) const;
    
protected:
    unsigned int seed;
    double standard_deviation;
    std::vector<double> displacement;
};

}
}
}

#endif
