/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SOLVER_SOLVER_HPP
#define SGPP_OPT_SLE_SOLVER_SOLVER_HPP

#include "opt/sle/system/System.hpp"

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
 * Abstract class for solving systems of linear equations.
 */
class Solver
{
public:
    /**
     * Constructor.
     */
    Solver()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~Solver()
    {
    }
    
    /**
     * Pure virtual method for a solving linear system.
     * 
     * @param       system  system to be solved
     * @param       b       right-hand side
     * @param[out]  x       solution to the system
     * @return              whether all went well (false if errors occurred)
     */
    virtual bool solve(system::System &system, const std::vector<double> &b,
                       std::vector<double> &x) const = 0;
    
    /**
     * Virtual method for solving multiple linear systems with different right-hand sides.
     * Defaults to calling the solve() method for a single right-hand side multiple times.
     * 
     * @param       system  system to be solved
     * @param       B       vector of right-hand sides
     * @param[out]  X       vector of solutions to the systems
     * @return              whether all went well (false if errors occurred)
     */
    virtual bool solve(system::System &system, const std::vector<std::vector<double> > &B,
                       std::vector<std::vector<double> > &X) const
    {
        std::vector<double> x;
        X.clear();
        
        for (size_t i = 0; i < B.size(); i++)
        {
            const std::vector<double> &b = B[i];
            
            if (solve(system, b, x))
            {
                X.push_back(x);
            } else
            {
                X.clear();
                return false;
            }
        }
        
        return true;
    }
};

}
}
}
}

#endif
