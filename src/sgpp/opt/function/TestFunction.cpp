#include "opt/function/TestFunction.hpp"

#include <random>
#include <iostream>

namespace sg
{
namespace opt
{
namespace function
{

const double TestFunction::DEFAULT_STANDARD_DEVIATION = 0.01;

TestFunction::TestFunction(size_t d) :
    ObjectiveFunction(d),
    displacement(std::vector<double>(d, 0.0))
{
}

double TestFunction::eval(const std::vector<double> &x)
{
    std::vector<double> x_displaced = x;
    displaceVector(x_displaced);
    return evalUndisplaced(x_displaced);
}

double TestFunction::getOptimalPoint(std::vector<double> &x)
{
    double fx = getOptimalPointUndisplaced(x);
    reverseDisplaceVector(x);
    return fx;
}

void TestFunction::generateDisplacement()
{
    generateDisplacement(std::random_device()());
}

void TestFunction::generateDisplacement(unsigned int seed)
{
    generateDisplacement(seed, DEFAULT_STANDARD_DEVIATION);
}

void TestFunction::generateDisplacement(unsigned int seed, double standard_deviation)
{
    this->seed = seed;
    this->standard_deviation = standard_deviation;
    
    std::mt19937 rng(seed);
    std::normal_distribution<> normal_dist(0, standard_deviation);
    
    displacement = std::vector<double>(d, 0.0);
    
    for (size_t t = 0; t < d; t++)
    {
        displacement[t] = normal_dist(rng);
    }
}

void TestFunction::displaceVector(std::vector<double> &x) const
{
    if (x.size() != d)
    {
        x.clear();
        std::cerr << "sg::opt::TestFunction::displaceVector: x doesn't match size\n";
    }
    
    for (size_t t = 0; t < d; t++)
    {
        x[t] += displacement[t];
    }
}

void TestFunction::reverseDisplaceVector(std::vector<double> &x) const
{
    if (x.size() != d)
    {
        x.clear();
        std::cerr << "sg::opt::TestFunction::reverseDisplaceVector: x doesn't match size\n";
    }
    
    for (size_t t = 0; t < d; t++)
    {
        x[t] -= displacement[t];
    }
}

unsigned int TestFunction::getSeed() const
{
    return seed;
}

double TestFunction::getStandardDeviation() const
{
    return standard_deviation;
}

void TestFunction::getDisplacement(std::vector<double> &displacement) const
{
    displacement = this->displacement;
}

}
}
}
