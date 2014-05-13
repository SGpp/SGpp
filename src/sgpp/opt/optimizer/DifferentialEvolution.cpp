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
const double DifferentialEvolution::DEFAULT_AVG_IMPROVEMENT_THRESHOLD = 1e-6;
const double DifferentialEvolution::DEFAULT_MAX_DISTANCE_THRESHOLD = 1e-4;

DifferentialEvolution::DifferentialEvolution(function::Objective &f) :
    DifferentialEvolution(f, DEFAULT_MAX_IT_COUNT)
{
}

DifferentialEvolution::DifferentialEvolution(function::Objective &f, size_t max_it_count) :
    DifferentialEvolution(f, max_it_count, std::random_device()())
{
}

DifferentialEvolution::DifferentialEvolution(function::Objective &f, size_t max_it_count,
                                             unsigned int seed) :
    DifferentialEvolution(f, max_it_count, seed, 10*f.getDimension(),
                          DEFAULT_CROSSOVER_PROBABILITY, DEFAULT_SCALING_FACTOR,
                          DEFAULT_IDLE_GENERATIONS_COUNT, DEFAULT_AVG_IMPROVEMENT_THRESHOLD,
                          DEFAULT_MAX_DISTANCE_THRESHOLD)
{
}

DifferentialEvolution::DifferentialEvolution(function::Objective &f, size_t max_it_count,
        unsigned int seed, size_t points_count,
        double crossover_probability, double scaling_factor,
        size_t idle_generations_count, double avg_improvement_threshold,
        double max_distance_threshold) :
    Optimizer(f, max_it_count),
    points_count(points_count),
    seed(seed),
    crossover_probability(crossover_probability),
    scaling_factor(scaling_factor),
    idle_generations_count(idle_generations_count),
    avg_improvement_threshold(avg_improvement_threshold),
    max_distance_threshold(max_distance_threshold)
{
}

void DifferentialEvolution::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (differential evolution)...");
    
    size_t d = f->getDimension();
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
        
        fx[i] = f->eval((*x_old)[i]);
    }
    
    double fopt = INFINITY;
    size_t xopt_index = 0;
    size_t last_nonidle_k = 0;
    double avg = 0.0;
    double last_avg = 0.0;
    size_t max_k = std::max(static_cast<size_t>(2), N / points_count) - 1;
    
    std::vector<std::vector<size_t> > a(max_k, std::vector<size_t>(points_count, 0)),
            b = a, c = a, j = a;
    std::vector<std::vector<std::vector<double> > > prob(max_k,
            std::vector<std::vector<double> >(points_count, std::vector<double>(d, 0)));
    
    for (size_t k = 0; k < max_k; k++)
    {
        for (size_t i = 0; i < points_count; i++)
        {
            do
            {
                a[k][i] = discr_dist_points(rng);
            } while (a[k][i] == i);
            
            do
            {
                b[k][i] = discr_dist_points(rng);
            } while ((b[k][i] == i) || (b[k][i] == a[k][i]));
            
            do
            {
                c[k][i] = discr_dist_points(rng);
            } while ((c[k][i] == i) || (c[k][i] == a[k][i]) || (c[k][i] == b[k][i]));
            
            j[k][i] = discr_dist_dim(rng);
            
            for (size_t t = 0; t < d; t++)
            {
                if (t != j[k][i])
                {
                    prob[k][i][t] = real_dist(rng);
                }
            }
        }
    }
    
    for (size_t k = 0; k < max_k; k++)
    {
        const std::vector<size_t> &a_k = a[k], &b_k = b[k], &c_k = c[k];
        const std::vector<size_t> &j_k = j[k];
        const std::vector<std::vector<double> > &prob_k = prob[k];
        
        #pragma omp parallel shared(k, a_k, b_k, c_k, j_k, prob_k, \
                                    x_old, d, fx, fopt, xopt_index, x_new) default(none)
        {
            std::vector<double> y(d, 0.0);
            std::unique_ptr<function::Objective> cur_f(f->clone());
            
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < points_count; i++)
            {
                const size_t &cur_a = a_k[i], &cur_b = b_k[i], &cur_c = c_k[i];
                const size_t &cur_j = j_k[i];
                const std::vector<double> &prob_ki = prob_k[i];
                bool in_domain = true;
                
                for (size_t t = 0; t < d; t++)
                {
                    const double &cur_prob = prob_ki[t];
                    
                    if ((t == cur_j) || (cur_prob < crossover_probability))
                    {
                        y[t] = (*x_old)[cur_a][t] +
                                scaling_factor * ((*x_old)[cur_b][t] - (*x_old)[cur_c][t]);
                    } else
                    {
                        y[t] = (*x_old)[i][t];
                    }
                    
                    if ((y[t] < 0.0) || (y[t] > 1.0))
                    {
                        in_domain = false;
                    }
                }
                
                double fy = (in_domain ? cur_f->eval(y) : INFINITY);
                
                if (fy < fx[i])
                {
                    #pragma omp critical
                    {
                        fx[i] = fy;
                        
                        if (fy < fopt)
                        {
                            xopt_index = i;
                            fopt = fy;
                        }
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
        }
        
        std::swap(x_old, x_new);
        avg = 0.0;
        
        for (size_t i = 0; i < points_count; i++)
        {
            avg += fx[i];
        }
        
        avg /= static_cast<double>(points_count);
        
        if (last_avg - avg >= avg_improvement_threshold)
        {
            last_nonidle_k = k;
        } else if (k - last_nonidle_k >= idle_generations_count)
        {
            double max_distance2 = 0.0;
            
            for (size_t i = 0; i < points_count; i++)
            {
                if (i == xopt_index)
                {
                    continue;
                }
                
                double distance2 = 0.0;
                
                for (size_t t = 0; t < d; t++)
                {
                    distance2 += ((*x_old)[i][t] - (*x_old)[xopt_index][t]) *
                                 ((*x_old)[i][t] - (*x_old)[xopt_index][t]);
                }
                
                if (distance2 > max_distance2)
                {
                    max_distance2 = distance2;
                }
            }
            
            if (std::sqrt(max_distance2) < max_distance_threshold)
            {
                break;
            }
        }
        
        last_avg = avg;
        
        if (k % 10 == 0)
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << fopt;
            tools::printer.printStatusUpdate(msg.str());
        }
    }
    
    xopt = (*x_old)[xopt_index];
    
    {
        std::stringstream msg;
        msg << max_k << " steps, f(x) = " << fopt;
        tools::printer.printStatusUpdate(msg.str());
    }
    
    tools::printer.printStatusEnd();
}

std::unique_ptr<Optimizer> DifferentialEvolution::clone()
{
    std::unique_ptr<Optimizer> result(new DifferentialEvolution(
            *f, N, seed, points_count, crossover_probability, scaling_factor,
            idle_generations_count, avg_improvement_threshold, max_distance_threshold));
    result->setStartingPoint(x0);
    return result;
}

}
}
}
