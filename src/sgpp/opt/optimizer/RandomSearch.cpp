#include "opt/optimizer/RandomSearch.hpp"
#include "opt/tools/Printer.hpp"

#include <algorithm>
#include <iostream>

namespace sg
{
namespace opt
{
namespace optimizer
{

RandomSearch::RandomSearch(function::ObjectiveFunction &f) :
    RandomSearch(f, DEFAULT_MAX_IT_COUNT)
{
}

RandomSearch::RandomSearch(function::ObjectiveFunction &f, size_t max_it_count) :
    RandomSearch(f, max_it_count, std::random_device()())
{
}

RandomSearch::RandomSearch(function::ObjectiveFunction &f, size_t max_it_count,
                           unsigned int seed) :
    RandomSearch(f, max_it_count, seed, default_optimizer,
                 std::min(10*f.getDimension(), static_cast<size_t>(100)))
{
}

RandomSearch::RandomSearch(function::ObjectiveFunction &f, size_t max_it_count, unsigned int seed,
                           Optimizer &optimizer, size_t points_count) :
    Optimizer(f, max_it_count),
    optimizer(optimizer),
    points_count(points_count),
    seed(seed),
    default_optimizer(NelderMead(f))
{
}

void RandomSearch::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (random search)...");
    
    size_t d = f.getDimension();
    std::vector<double> x0(d, 0.0);
    
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    size_t remaining_N = N;
    std::vector<double> cur_xopt(d, 0.0);
    
    double cur_fopt = INFINITY;
    double fopt = INFINITY;
    
    for (size_t k = 0; k < points_count; k++)
    {
        for (size_t t = 0; t < d; t++)
        {
            x0[t] = dist(rng);
        }
        
        size_t cur_N = static_cast<size_t>(std::ceil(static_cast<double>(remaining_N) /
                                             static_cast<double>(points_count - k)));
        
        optimizer.setStartingPoint(x0);
        optimizer.setMaxItCount(cur_N);
        
        tools::printer.increaseCurrentLevel(2);
        optimizer.optimize(cur_xopt);
        tools::printer.decreaseCurrentLevel(2);
        
        remaining_N -= cur_N;
        cur_fopt = f.eval(cur_xopt);
        
        if (cur_fopt < fopt)
        {
            fopt = cur_fopt;
            xopt = cur_xopt;
        }
        
        {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(k) / static_cast<double>(points_count) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             ", f(x) = " + std::to_string(fopt));
        }
    }
    
    {
        std::stringstream msg;
        msg << "100.0%, f(x) = " << fopt;
        tools::printer.printStatusUpdate(msg.str());
    }
    
    tools::printer.printStatusEnd();
}

}
}
}
