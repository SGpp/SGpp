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
    static const size_t DEFAULT_MAX_IT_COUNT = 1000;
    
    RandomSearch(function::Objective &f, size_t max_it_count = DEFAULT_MAX_IT_COUNT);
    RandomSearch(function::Objective &f, size_t max_it_count, unsigned int seed);
    RandomSearch(size_t max_it_count, unsigned int seed,
                 Optimizer &optimizer, size_t points_count);
    
    std::unique_ptr<Optimizer> clone();
    
    double optimize(std::vector<double> &xopt);
    
protected:
    std::unique_ptr<Optimizer> optimizer;
    Optimizer *optimizer_to_use;
    size_t points_count;
    unsigned int seed;
    
    NelderMead default_optimizer;
    
    void initialize(unsigned int seed, size_t points_count);
};

}
}
}

#endif
