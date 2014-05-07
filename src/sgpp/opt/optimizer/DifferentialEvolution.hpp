#ifndef SGPP_OPT_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP
#define SGPP_OPT_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP

#include "opt/optimizer/Optimizer.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

class DifferentialEvolution : public Optimizer
{
public:
    static const double DEFAULT_CROSSOVER_PROBABILITY;
    static const double DEFAULT_SCALING_FACTOR;
    
    DifferentialEvolution(function::ObjectiveFunction &f);
    DifferentialEvolution(function::ObjectiveFunction &f, size_t max_it_count);
    DifferentialEvolution(function::ObjectiveFunction &f, size_t max_it_count, unsigned int seed);
    DifferentialEvolution(function::ObjectiveFunction &f, size_t max_it_count,
                          unsigned int seed, size_t points_count,
                          double crossover_probability, double scaling_factor);
    
    void optimize(std::vector<double> &xopt);
    
protected:
    size_t points_count;
    size_t seed;
    double crossover_probability;
    double scaling_factor;
};

}
}
}

#endif
