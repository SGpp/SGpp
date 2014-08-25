/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/optimizer/NelderMead.hpp"
#include "opt/tools/Permuter.hpp"
#include "opt/tools/Printer.hpp"

#include <algorithm>
#include <iostream>

namespace sg
{
namespace opt
{
namespace optimizer
{

const double NelderMead::DEFAULT_ALPHA = 1.0;
const double NelderMead::DEFAULT_BETA = 2.0;
const double NelderMead::DEFAULT_GAMMA = 0.5;
const double NelderMead::DEFAULT_DELTA = 0.5;
const double NelderMead::STARTING_SIMPLEX_EDGE_LENGTH = 0.4;

NelderMead::NelderMead(function::Objective &f,
        size_t max_fcn_eval_count, double alpha, double beta, double gamma, double delta) :
    Optimizer(f, max_fcn_eval_count),
    alpha(alpha),
    beta(beta),
    gamma(gamma),
    delta(delta)
{
}

double NelderMead::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (Nelder-Mead)...");
    
    size_t d = f->getDimension();
    std::vector<std::vector<double> > points(d+1, x0);
    std::vector<std::vector<double> > points_new(d+1, x0);
    std::vector<double> f_points(d+1, 0.0);
    std::vector<double> f_points_new(d+1, 0.0);
    
    // construct starting simplex
    for (size_t t = 0; t < d; t++)
    {
        points[t+1][t] = std::min(points[t+1][t] + STARTING_SIMPLEX_EDGE_LENGTH, 1.0);
        f_points[t+1] = f->eval(points[t+1]);
    }
    
    f_points[0] = f->eval(points[0]);
    
    std::vector<size_t> index(d+1, 0);
    std::vector<double> point_o(d, 0.0);
    std::vector<double> point_r(d, 0.0);
    std::vector<double> point_e(d, 0.0);
    std::vector<double> point_ic(d, 0.0);
    std::vector<double> point_oc(d, 0.0);
    size_t k = 0;
    size_t number_of_fcn_evals = d+1;
    
    while (true)
    {
        // sort points by function value
        for (size_t i = 0; i < d+1; i++) {
            index[i] = i;
        }
        
        {
            tools::Permuter<double> permuter(f_points);
            std::sort(index.begin(), index.end(), permuter);
        }
        
        // that could be solved more efficiently, but it suffices for now
        for (size_t i = 0; i < d+1; i++) {
            points_new[i] = points[index[i]];
            f_points_new[i] = f_points[index[i]];
        }
        
        points = points_new;
        f_points = f_points_new;
        
        bool in_domain = true;
        bool shrink = false;
        
        // calculate point_o (barycenter of all points but the last) and
        // point_r (reflected point) simultaneously
        for (size_t t = 0; t < d; t++)
        {
            point_o[t] = 0.0;
            
            for (size_t i = 0; i < d; i++)
            {
                point_o[t] += points[i][t];
            }
            
            point_o[t] /= static_cast<double>(d);
            point_r[t] = point_o[t] + alpha * (point_o[t] - points[d][t]);
            
            if ((point_r[t] < 0.0) || (point_r[t] > 1.0))
            {
                in_domain = false;
            }
        }
        
        double f_point_r = (in_domain ? f->eval(point_r) : INFINITY);
        number_of_fcn_evals++;
        
        if ((f_points[0] <= f_point_r) && (f_point_r < f_points[d-1]))
        {
            points[d] = point_r;
            f_points[d] = f_point_r;
        } else if (f_point_r < f_points[0])
        {
            bool in_domain = true;
            
            // calculate expanded point
            for (size_t t = 0; t < d; t++)
            {
                point_e[t] = point_o[t] + beta * (point_r[t] - point_o[t]);
                
                if ((point_e[t] < 0.0) || (point_e[t] > 1.0))
                {
                    in_domain = false;
                }
            }
            
            double f_point_e = (in_domain ? f->eval(point_e) : INFINITY);
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
        } else if (f_point_r < f_points[d])
        {
            bool in_domain = true;
            
            // calculate outer contracted point
            for (size_t t = 0; t < d; t++)
            {
                point_oc[t] = point_o[t] + gamma * (point_r[t] - point_o[t]);
                
                if ((point_oc[t] < 0.0) || (point_oc[t] > 1.0))
                {
                    in_domain = false;
                }
            }
            
            double f_point_oc = (in_domain ? f->eval(point_oc) : INFINITY);
            number_of_fcn_evals++;
            
            if (f_point_oc <= f_point_r)
            {
                points[d] = point_oc;
                f_points[d] = f_point_oc;
            } else
            {
                shrink = true;
            }
        } else
        {
            bool in_domain = true;
            
            // calculate inner contracted point
            for (size_t t = 0; t < d; t++)
            {
                point_ic[t] = point_o[t] - gamma * (point_o[t] - points[d][t]);
                
                if ((point_ic[t] < 0.0) || (point_ic[t] > 1.0))
                {
                    in_domain = false;
                }
            }
            
            double f_point_ic = (in_domain ? f->eval(point_ic) : INFINITY);
            number_of_fcn_evals++;
            
            if (f_point_ic < f_points[d])
            {
                points[d] = point_ic;
                f_points[d] = f_point_ic;
            } else
            {
                shrink = true;
            }
        }
        
        if (shrink)
        {
            // shrink all points but the first
            for (size_t i = 1; i < d+1; i++)
            {
                bool in_domain = true;
                
                for (size_t t = 0; t < d; t++)
                {
                    points[i][t] = points[0][t] + delta * (points[i][t] - points[0][t]);
                    
                    if ((points[i][t] < 0.0) || (points[i][t] > 1.0))
                    {
                        in_domain = false;
                    }
                }
                
                f_points[i] = (in_domain ? f->eval(points[i]) : INFINITY);
            }
            
            number_of_fcn_evals += d;
        }
        
        // status printing
        if (k % 10 == 0)
        {
            tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " +
                                             toString(f_points[0]));
        }
        
        if (number_of_fcn_evals + (d+2) > N)
        {
            break;
        }
        
        k++;
    }
    
    xopt = points[0];
    
    tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(f_points[0]));
    tools::printer.printStatusEnd();
    
    return f_points[0];
}

tools::SmartPointer<Optimizer> NelderMead::clone()
{
    tools::SmartPointer<Optimizer> result(new NelderMead(*f, N, alpha, beta, gamma, delta));
    result->setStartingPoint(x0);
    return result;
}

double NelderMead::getAlpha() const
{
    return alpha;
}

void NelderMead::setAlpha(double alpha)
{
    this->alpha = alpha;
}

double NelderMead::getBeta() const
{
    return beta;
}

void NelderMead::setBeta(double beta)
{
    this->beta = beta;
}

double NelderMead::getGamma() const
{
    return gamma;
}

void NelderMead::setGamma(double gamma)
{
    this->gamma = gamma;
}

double NelderMead::getDelta() const
{
    return delta;
}

void NelderMead::setDelta(double delta)
{
    this->delta = delta;
}

}
}
}
