/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGSOLVERINFO_HPP
#define SGSOLVERINFO_HPP

namespace sg
{
/**
 * This struct is needed for exchanging the SGSolver Settings
 * to another address space.
 */
struct SGSolverInfo
{
	/// lambda, the regression parameter
	double lambda;
	/// final error of SGSovler
	double epsilon;
	/// maximum number of iterations
	size_t imax;
	/// StiffnessMode
	std::string StiffnessMode;
	/// GridType
	std::string GridType;
};

}

#endif /* SGSOLVERINFO_HPP */
