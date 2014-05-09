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
    static const size_t DEFAULT_IDLE_GENERATIONS_COUNT = 20;
    static const double DEFAULT_AVG_IMPROVEMENT_THRESHOLD;
    static const double DEFAULT_MAX_DISTANCE_THRESHOLD;
    
    DifferentialEvolution(function::Objective &f);
    DifferentialEvolution(function::Objective &f, size_t max_it_count);
    DifferentialEvolution(function::Objective &f, size_t max_it_count, unsigned int seed);
    DifferentialEvolution(function::Objective &f, size_t max_it_count,
                          unsigned int seed, size_t points_count,
                          double crossover_probability, double scaling_factor,
                          size_t idle_generations_count,
                          double avg_improvement_threshold,
                          double max_distance_threshold);
    
    void optimize(std::vector<double> &xopt);
    
protected:
    size_t points_count;
    size_t seed;
    double crossover_probability;
    double scaling_factor;
    size_t idle_generations_count;
    double avg_improvement_threshold;
    double max_distance_threshold;
};

}
}
}

#endif
