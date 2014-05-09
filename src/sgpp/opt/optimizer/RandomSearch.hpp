#ifndef SGPP_OPT_OPTIMIZER_RANDOMSEARCH_HPP
#define SGPP_OPT_OPTIMIZER_RANDOMSEARCH_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/optimizer/NelderMead.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

class RandomSearch : public Optimizer
{
public:
    RandomSearch(function::Objective &f);
    RandomSearch(function::Objective &f, size_t max_it_count);
    RandomSearch(function::Objective &f, size_t max_it_count, unsigned int seed);
    RandomSearch(function::Objective &f, size_t max_it_count, unsigned int seed,
                 Optimizer &optimizer, size_t points_count);
    
    void optimize(std::vector<double> &xopt);
    
protected:
    Optimizer &optimizer;
    size_t points_count;
    size_t seed;
    
    NelderMead default_optimizer;
};

}
}
}

#endif
