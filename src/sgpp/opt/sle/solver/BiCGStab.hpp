/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SOLVER_BICGSTAB_HPP
#define SGPP_OPT_SLE_SOLVER_BICGSTAB_HPP

#include "opt/sle/solver/Solver.hpp"

#include <cstddef>
#include <vector>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

/**
 * Linear system solver implementing the iterative BiCGStab method.
 */
class BiCGStab : public Solver
{
public:
    /// default maximal number of iterations
    static const size_t DEFAULT_MAX_IT_COUNT = 1000;
    /// default tolerance
    static const double DEFAULT_TOLERANCE;
    
    /**
     * Constructor.
     */
    BiCGStab();
    
    /**
     * @param max_it_count      maximal number of iterations
     * @param tolerance         tolerance
     * @param starting_point    starting vector
     */
    BiCGStab(size_t max_it_count, double tolerance, const std::vector<double> &starting_point);
    
    /**
     * @param       system  system to be solved
     * @param       b       right-hand side
     * @param[out]  x       solution to the system
     * @return              whether all went well (false if errors occurred)
     */
    bool solve(system::System &system, const std::vector<double> &b, std::vector<double> &x) const;
    
    /**
     * @return              maximal number of iterations
     */
    size_t getMaxItCount() const;
    
    /**
     * @param max_it_count  maximal number of iterations
     */
    void setMaxItCount(size_t max_it_count);
    
    /**
     * @return              tolerance
     */
    double getTolerance() const;
    
    /**
     * @param tolerance     tolerance
     */
    void setTolerance(double tolerance);
    
    /**
     * @return                  starting vector
     */
    const std::vector<double> &getStartingPoint() const;
    
    /**
     * @param starting_point    starting vector
     */
    void setStartingPoint(const std::vector<double> &starting_point);
    
protected:
    /// maximal number of iterations
    size_t N;
    /// tolerance
    double tol;
    /// starting vector
    std::vector<double> x0;
};

}
}
}
}

#endif
