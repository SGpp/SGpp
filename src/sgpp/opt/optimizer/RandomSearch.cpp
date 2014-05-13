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

RandomSearch::RandomSearch(function::Objective &f) :
    RandomSearch(f, DEFAULT_MAX_IT_COUNT)
{
}

RandomSearch::RandomSearch(function::Objective &f, size_t max_it_count) :
    RandomSearch(f, max_it_count, std::random_device()())
{
}

RandomSearch::RandomSearch(function::Objective &f, size_t max_it_count,
                           unsigned int seed) :
    Optimizer(f, max_it_count),
    optimizer(nullptr),
    optimizer_to_use(&default_optimizer),
    points_count(std::min(10*f.getDimension(), static_cast<size_t>(100))),
    seed(seed),
    default_optimizer(NelderMead(f))
{
}

RandomSearch::RandomSearch(size_t max_it_count, unsigned int seed,
                           Optimizer &optimizer, size_t points_count) :
    Optimizer(*optimizer.getObjectiveFunction(), max_it_count),
    optimizer(optimizer.clone()),
    optimizer_to_use(this->optimizer.get()),
    points_count(points_count),
    seed(seed),
    default_optimizer(NelderMead(*f))
{
}

void RandomSearch::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (random search)...");
    
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    size_t d = f->getDimension();
    std::vector<std::vector<double> > x0(points_count, std::vector<double>(d, 0.0));
    std::vector<size_t> round_N(points_count, 0);
    size_t remaining_N = N;
    
    for (size_t k = 0; k < points_count; k++)
    {
        round_N[k] = static_cast<size_t>(std::ceil(static_cast<double>(remaining_N) /
                                         static_cast<double>(points_count - k)));
        remaining_N -= round_N[k];
        
        for (size_t t = 0; t < d; t++)
        {
            x0[k][t] = dist(rng);
        }
    }
    
    double fopt = INFINITY;
    
    tools::printer.disableStatusPrinting();
    
    #pragma omp parallel shared(d, x0, round_N, xopt, fopt, tools::printer) default(none)
    {
        std::unique_ptr<Optimizer> cur_optimizer(optimizer_to_use->clone());
        
        std::vector<double> cur_xopt(d, 0.0);
        double cur_fopt;
        
        #pragma omp for ordered schedule(dynamic)
        for (size_t k = 0; k < points_count; k++)
        {
            //optimizer.setObjectiveFunction(f);
            cur_optimizer->setStartingPoint(x0[k]);
            cur_optimizer->setMaxItCount(round_N[k]);
            cur_optimizer->optimize(cur_xopt);
            
            cur_fopt = cur_optimizer->getObjectiveFunction()->eval(cur_xopt);
            
            #pragma omp critical
            {
                if (cur_fopt < fopt)
                {
                    fopt = cur_fopt;
                    xopt = cur_xopt;
                }
            }
            
            #pragma omp ordered
            {
                char str[10];
                snprintf(str, 10, "%.1f%%",
                         static_cast<double>(k) / static_cast<double>(points_count) * 100.0);
                tools::printer.getMutex().lock();
                tools::printer.enableStatusPrinting();
                tools::printer.printStatusUpdate(std::string(str) +
                                                 ", f(x) = " + std::to_string(fopt));
                tools::printer.disableStatusPrinting();
                tools::printer.getMutex().unlock();
            }
        }
    }
    
    tools::printer.enableStatusPrinting();
    
    {
        std::stringstream msg;
        msg << "100.0%, f(x) = " << fopt;
        tools::printer.printStatusUpdate(msg.str());
    }
    
    tools::printer.printStatusEnd();
}

std::unique_ptr<Optimizer> RandomSearch::clone()
{
    std::unique_ptr<Optimizer> result(new RandomSearch(N, seed, *optimizer_to_use, points_count));
    result->setStartingPoint(x0);
    return result;
}

}
}
}
