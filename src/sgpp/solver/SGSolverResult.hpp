/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGSOLVERRESULT_HPP
#define SGSOLVERRESULT_HPP

namespace sg
{
namespace solver
{
/**
 * This struct is needed for exchanging the SGSolver Results
 * to another address space.
 */
struct SGSolverResult
{
	/// Iterations that where needed to solve the system
	size_t nNeededIterations;
	/// Final residuum of SGSolver
	double finalResiduum;
};

}
}

#endif /* SGSOLVERRESULT_HPP */
