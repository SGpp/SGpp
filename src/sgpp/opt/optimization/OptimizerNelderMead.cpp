#include "opt/optimization/OptimizerNelderMead.hpp"
#include "opt/tools/Printer.hpp"

#include <algorithm>
#include <iostream>

namespace sg
{
namespace opt
{
namespace optimization
{
    
const double OptimizerNelderMead::DEFAULT_ALPHA = 1.0;
const double OptimizerNelderMead::DEFAULT_GAMMA = 2.0;
const double OptimizerNelderMead::DEFAULT_RHO = -0.5;
const double OptimizerNelderMead::DEFAULT_SIGMA = 0.5;
const double OptimizerNelderMead::STARTING_SIMPLEX_EDGE_LENGTH = 0.4;

OptimizerNelderMead::OptimizerNelderMead(function::ObjectiveFunction &f) :
    OptimizerNelderMead(f, DEFAULT_MAX_IT_COUNT,
                        DEFAULT_ALPHA, DEFAULT_GAMMA, DEFAULT_RHO, DEFAULT_SIGMA,
                        DEFAULT_MAX_FCN_EVAL_COUNT)
{
}

OptimizerNelderMead::OptimizerNelderMead(
        function::ObjectiveFunction &f,
        size_t max_it_count, double alpha, double gamma, double rho, double sigma,
        size_t max_fcn_eval_count) :
    Optimizer(f, max_it_count),
    alpha(alpha),
    gamma(gamma),
    rho(rho),
    sigma(sigma),
    max_fcn_eval_count(max_fcn_eval_count)
{
}

void OptimizerNelderMead::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (Nelder-Mead)...");
    
    size_t d = f.getDimension();
    std::vector<std::vector<double> > points(d+1, x0);
    std::vector<std::vector<double> > points_new(d+1, x0);
    std::vector<double> f_points(d+1, 0.0);
    std::vector<double> f_points_new(d+1, 0.0);
    
    for (size_t t = 0; t < d; t++)
    {
        //points[0][t] = std::max(points[0][t] - 0.1, 0.0);
        //points[t+1][t] = std::min(points[t+1][t] + 0.1, 1.0);
        points[t+1][t] = std::min(points[t+1][t] + STARTING_SIMPLEX_EDGE_LENGTH, 1.0);
        f_points[t+1] = f.eval(points[t+1]);
        //std::cout << "i = " << i+1 << ": " << points[i+1] << ", " << f_points[i+1] << ", " << f(points[i+1], data) << "\n";
    }
    
    f_points[0] = f.eval(points[0]);
    
    std::vector<size_t> index(d+1, 0);
    std::vector<double> point_0(d, 0.0);
    std::vector<double> point_r(d, 0.0);
    std::vector<double> point_e(d, 0.0);
    std::vector<double> point_c(d, 0.0);
    size_t k = 0;
    size_t number_of_fcn_evals = d+1;
    
    for (k = 0; k < N; k++)
    {
        /*std::cout << "k = " << k << ", number_of_fcn_evals = " << number_of_fcn_evals << "\n";
        for (size_t i = 0; i < d+1; i++)
        {
            std::cout << "points[" << i << "] = [" << points[i][0] << ", " << points[i][1] << "]\n";
            std::cout << "f_points[" << i << "] = [" << f_points[i] << "]\n";
        }*/
        
        for (size_t i = 0; i < d+1; i++) {
            index[i] = i;
        }
        
        std::sort(index.begin(), index.end(),
                [&](const size_t &a, const size_t &b) {
                    return (f_points[a] < f_points[b]);
                }
        );
        
        for (size_t i = 0; i < d+1; i++) {
            points_new[i] = points[index[i]];
            f_points_new[i] = f_points[index[i]];
            //std::cout << "i = " << i << ": " << points_new[i] << ", " << f_points_new[i] << ", " << f(points_new[i], data) << "\n";
        }
        
        points = points_new;
        f_points = f_points_new;
        
        //bool converged = true;
        bool in_domain = true;
        
        for (size_t t = 0; t < d; t++)
        {
            point_0[t] = 0.0;
            //double point_min = 1.0;
            //double point_max = 0.0;
            
            for (size_t i = 0; i < d; i++)
            {
                point_0[t] += points[i][t];
                
                /*if (points[j][i] < point_min)
                {
                    point_min = points[j][i];
                }
                
                if (points[j][i] > point_max)
                {
                    point_max = points[j][i];
                }*/
            }
            
            /*if (point_max - point_min >= tol)
            {
                converged = false;
            }*/
            
            point_0[t] /= (double)d;
            point_r[t] = point_0[t] + alpha * (point_0[t] - points[d][t]);
            
            if ((point_r[t] < 0.0) || (point_r[t] > 1.0))
            {
                in_domain = false;
            }
        }
        
        /*if (converged)
        {
            break;
        }*/
        
        double f_point_r = (in_domain ? f.eval(point_r) : INFINITY);
        number_of_fcn_evals++;
        
        if ((f_points[0] <= f_point_r) && (f_point_r < f_points[d-1]))
        {
            points[d] = point_r;
            f_points[d] = f_point_r;
        } else if (f_point_r < f_points[0])
        {
            in_domain = true;
            
            for (size_t t = 0; t < d; t++)
            {
                point_e[t] = point_0[t] + gamma * (point_0[t] - points[d][t]);
                
                if ((point_e[t] < 0.0) || (point_e[t] > 1.0))
                {
                    in_domain = false;
                }
            }
            
            double f_point_e = (in_domain ? f.eval(point_e) : INFINITY);
            number_of_fcn_evals++;
            
            if (f_point_e < f_point_r)
            {
                points[d] = point_e;
                f_points[d] = f_point_e;
            } else
            {
                points[d] = point_r;
                f_points[d] = f_point_r;
            }
        } else
        {
            in_domain = true;
            
            for (size_t t = 0; t < d; t++)
            {
                point_c[t] = point_0[t] + rho * (point_0[t] - points[d][t]);
                
                if ((point_c[t] < 0.0) || (point_c[t] > 1.0))
                {
                    in_domain = false;
                }
            }
            
            double f_point_c = (in_domain ? f.eval(point_c) : INFINITY);
            number_of_fcn_evals++;
            
            if (f_point_c < f_points[d])
            {
                points[d] = point_c;
                f_points[d] = f_point_c;
            } else
            {
                for (size_t i = 1; i < d+1; i++)
                {
                    in_domain = true;
                    
                    for (size_t t = 0; t < d; t++)
                    {
                        points[i][t] = points[0][t] + sigma * (points[i][t] - points[0][t]);
                        
                        if ((points[i][t] < 0.0) || (points[i][t] > 1.0))
                        {
                            in_domain = false;
                        }
                    }
                    
                    f_points[i] = (in_domain ? f.eval(points[i]) : INFINITY);
                }
                
                number_of_fcn_evals += d;
            }
        }
        
        if (k % 10 == 0)
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << f_points[0];
            tools::printer.printStatusUpdate(msg.str());
        }
        
        if (number_of_fcn_evals + (d+2) > max_fcn_eval_count)
        {
            break;
        }
    }
    
    xopt = points[0];
    
    std::stringstream msg;
    msg << k << " steps, f(x) = " << f_points[0];
    tools::printer.printStatusUpdate(msg.str());
    tools::printer.printStatusEnd();
}

double OptimizerNelderMead::getAlpha() const
{
    return alpha;
}

void OptimizerNelderMead::setAlpha(double alpha)
{
    this->alpha = alpha;
}

double OptimizerNelderMead::getGamma() const
{
    return gamma;
}

void OptimizerNelderMead::setGamma(double gamma)
{
    this->gamma = gamma;
}

double OptimizerNelderMead::getRho() const
{
    return rho;
}

void OptimizerNelderMead::setRho(double rho)
{
    this->rho = rho;
}

double OptimizerNelderMead::getSigma() const
{
    return rho;
}

void OptimizerNelderMead::setSigma(double sigma)
{
    this->sigma = sigma;
}

size_t OptimizerNelderMead::getMaxFcnEvalCount() const
{
    return max_fcn_eval_count;
}

void OptimizerNelderMead::setMaxFcnEvalCount(size_t max_fcn_eval_count)
{
    this->max_fcn_eval_count = max_fcn_eval_count;
}

}
}
}
