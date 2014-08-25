/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/sle/solver/BiCGStab.hpp"
#include "opt/tools/Printer.hpp"

#include <cmath>
#include <numeric>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

const double BiCGStab::DEFAULT_TOLERANCE = 1e-10;

BiCGStab::BiCGStab() :
    Solver(),
    N(DEFAULT_MAX_IT_COUNT),
    tol(DEFAULT_TOLERANCE),
    x0(std::vector<double>())
{
}

BiCGStab::BiCGStab(size_t max_it_count, double tolerance, const std::vector<double> &x0) :
    Solver(),
    N(max_it_count),
    tol(tolerance),
    x0(x0)
{
}

bool BiCGStab::solve(system::System &system, const std::vector<double> &b,
                     std::vector<double> &x) const
{
    tools::printer.printStatusBegin("Solving linear system (BiCGStab)...");
    
    const size_t n = b.size();
    std::vector<double> r(n, 0.0);
    
    std::vector<double> my_x0 = x0;
    
    if (my_x0.empty())
    {
        my_x0 = std::vector<double>(n, 0.0);
    }
    
    x = std::vector<double>(n, 0.0);
    system.matrixVectorMultiplication(my_x0, r);
    
    for (size_t i = 0; i < n; i++)
    {
        r[i] = b[i] - r[i];
    }
    
    std::vector<double> r0hat(r);
    double rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    std::vector<double> v(n, 0.0);
    std::vector<double> p(n, 0.0);
    std::vector<double> s(n, 0.0);
    std::vector<double> t(n, 0.0);
    double r_norm_squared = 0.0;
    
    for (size_t k = 0; k < N; k++)
    {
        double last_rho = rho;
        rho = std::inner_product(r0hat.begin(), r0hat.end(), r.begin(), 0.0);
        double beta = (rho / last_rho) * (alpha / omega);
        
        for (size_t i = 0; i < n; i++)
        {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        
        system.matrixVectorMultiplication(p, v);
        alpha = rho / std::inner_product(r0hat.begin(), r0hat.end(), v.begin(), 0.0);
        
        for (size_t i = 0; i < n; i++)
        {
            s[i] = r[i] - alpha * v[i];
        }
        
        system.matrixVectorMultiplication(s, t);
        omega = std::inner_product(t.begin(), t.end(), s.begin(), 0.0) /
                std::inner_product(t.begin(), t.end(), t.begin(), 0.0);
        
        if (std::isnan(omega))
        {
            tools::printer.printStatusEnd("error: could not solve linear system!");
            return false;
        }
        
        for (size_t i = 0; i < n; i++)
        {
            x[i] = x[i] + alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }
        
        r_norm_squared = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
        
        tools::printer.printStatusUpdate("k = " + toString(k) + ", residual norm = " +
                                         toString(sqrt(r_norm_squared)));
        
        if (r_norm_squared < tol*tol)
        {
            break;
        }
    }
    
    tools::printer.printStatusUpdate("residual norm = " + toString(sqrt(r_norm_squared)));
    tools::printer.printStatusEnd();
    return true;
}

size_t BiCGStab::getMaxItCount() const
{
    return N;
}

void BiCGStab::setMaxItCount(size_t max_it_count)
{
    N = max_it_count;
}

double BiCGStab::getTolerance() const
{
    return tol;
}

void BiCGStab::setTolerance(double tolerance)
{
    tol = tolerance;
}

const std::vector<double> &BiCGStab::getStartingPoint() const
{
    return x0;
}

void BiCGStab::setStartingPoint(const std::vector<double> &starting_point)
{
    x0 = starting_point;
}

}
}
}
}
