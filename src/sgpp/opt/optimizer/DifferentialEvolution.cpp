#include "opt/optimizer/DifferentialEvolution.hpp"
#include "opt/tools/Printer.hpp"

#include <algorithm>
#include <iostream>

namespace sg
{
namespace opt
{
namespace optimizer
{

const double DifferentialEvolution::DEFAULT_CROSSOVER_PROBABILITY = 0.5;
const double DifferentialEvolution::DEFAULT_SCALING_FACTOR = 0.6;

DifferentialEvolution::DifferentialEvolution(function::ObjectiveFunction &f) :
    DifferentialEvolution(f, DEFAULT_MAX_IT_COUNT)
{
}

DifferentialEvolution::DifferentialEvolution(function::ObjectiveFunction &f, size_t max_it_count) :
    DifferentialEvolution(f, max_it_count, std::random_device()())
{
}

DifferentialEvolution::DifferentialEvolution(function::ObjectiveFunction &f, size_t max_it_count,
                                             unsigned int seed) :
    DifferentialEvolution(f, max_it_count, seed, 10*f.getDimension(),
                          DEFAULT_CROSSOVER_PROBABILITY, DEFAULT_SCALING_FACTOR)
{
}

DifferentialEvolution::DifferentialEvolution(function::ObjectiveFunction &f, size_t max_it_count,
        unsigned int seed, size_t points_count,
        double crossover_probability, double scaling_factor) :
    Optimizer(f, max_it_count),
    points_count(points_count),
    seed(seed),
    crossover_probability(crossover_probability),
    scaling_factor(scaling_factor)
{
}

void DifferentialEvolution::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (differential evolution)...");
    
    size_t d = f.getDimension();
    std::vector<std::vector<double> > x1(points_count, std::vector<double>(d, 0.0));
    std::vector<double> fx(points_count, 0.0);
    
    std::vector<std::vector<double> > x2 = x1;
    
    std::vector<std::vector<double> > *x_old = &x1;
    std::vector<std::vector<double> > *x_new = &x2;
    
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::uniform_int_distribution<size_t> discr_dist_points(0, points_count-1);
    std::uniform_int_distribution<size_t> discr_dist_dim(0, d-1);
    
    for (size_t i = 0; i < points_count; i++)
    {
        for (size_t t = 0; t < d; t++)
        {
            (*x_old)[i][t] = real_dist(rng);
        }
        
        fx[i] = f.eval((*x_old)[i]);
    }
    
    double fopt = INFINITY;
    size_t xopt_index = 0;
    std::vector<double> y(d, 0.0);
    size_t k = 0;
    size_t number_of_fcn_evals = points_count;
    
    while (true)
    {
        for (size_t i = 0; i < points_count; i++)
        {
            size_t a, b, c;
            bool in_domain = true;
            
            do
            {
                a = discr_dist_points(rng);
            } while (a == i);
            
            do
            {
                b = discr_dist_points(rng);
            } while ((b == i) || (b == a));
            
            do
            {
                c = discr_dist_points(rng);
            } while ((c == i) || (c == a) || (c == b));
            
            size_t j = discr_dist_dim(rng);
            
            for (size_t t = 0; t < d; t++)
            {
                if ((t == j) || (real_dist(rng) < crossover_probability))
                {
                    y[t] = (*x_old)[a][t] + scaling_factor * ((*x_old)[b][t] - (*x_old)[c][t]);
                } else
                {
                    y[t] = (*x_old)[i][t];
                }
                
                if ((y[t] < 0.0) || (y[t] > 1.0))
                {
                    in_domain = false;
                }
            }
            
            double fy = (in_domain ? f.eval(y) : INFINITY);
            
            if (fy < fx[i])
            {
                fx[i] = fy;
                
                if (fy < fopt)
                {
                    xopt_index = i;
                    fopt = fy;
                }
                
                for (size_t t = 0; t < d; t++)
                {
                    (*x_new)[i][t] = y[t];
                }
            } else
            {
                for (size_t t = 0; t < d; t++)
                {
                    (*x_new)[i][t] = (*x_old)[i][t];
                }
            }
        }
        
        std::swap(x_old, x_new);
        
        number_of_fcn_evals += points_count;
        
        if (number_of_fcn_evals + points_count > N)
        {
            break;
        }
        
        if (k % 10 == 0)
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << fopt;
            tools::printer.printStatusUpdate(msg.str());
        }
        
        k++;
    }
    
    xopt = (*x_old)[xopt_index];
    
    {
        std::stringstream msg;
        msg << k << " steps, f(x) = " << fopt;
        tools::printer.printStatusUpdate(msg.str());
    }
    
    tools::printer.printStatusEnd();
}

}
}
}
