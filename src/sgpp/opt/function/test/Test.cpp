#include "opt/function/test/Test.hpp"

#include <random>
#include <iostream>

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

const double Test::DEFAULT_STANDARD_DEVIATION = 0.01;

Test::Test(size_t d) :
    Objective(d),
    displacement(std::vector<double>(d, 0.0))
{
}

double Test::eval(const std::vector<double> &x)
{
    std::vector<double> x_displaced = x;
    displaceVector(x_displaced);
    return evalUndisplaced(x_displaced);
}

double Test::getOptimalPoint(std::vector<double> &x)
{
    double fx = getOptimalPointUndisplaced(x);
    reverseDisplaceVector(x);
    return fx;
}

void Test::generateDisplacement()
{
    generateDisplacement(std::random_device()());
}

void Test::generateDisplacement(unsigned int seed)
{
    generateDisplacement(seed, DEFAULT_STANDARD_DEVIATION);
}

void Test::generateDisplacement(unsigned int seed, double standard_deviation)
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

void Test::displaceVector(std::vector<double> &x) const
{
    if (x.size() != d)
    {
        x.clear();
        std::cerr << "sg::opt::Test::displaceVector: x doesn't match size\n";
        return;
    }
    
    for (size_t t = 0; t < d; t++)
    {
        x[t] += displacement[t];
    }
}

void Test::reverseDisplaceVector(std::vector<double> &x) const
{
    if (x.size() != d)
    {
        x.clear();
        std::cerr << "sg::opt::Test::reverseDisplaceVector: x doesn't match size\n";
        return;
    }
    
    for (size_t t = 0; t < d; t++)
    {
        x[t] -= displacement[t];
    }
}

unsigned int Test::getSeed() const
{
    return seed;
}

double Test::getStandardDeviation() const
{
    return standard_deviation;
}

void Test::getDisplacement(std::vector<double> &displacement) const
{
    displacement = this->displacement;
}

}
}
}
}
