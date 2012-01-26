/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LEARNER_HPP
#define LEARNER_HPP

namespace sg
{

namespace base
{

struct RegularGridConfiguration
{
	GirdType type_;
	size_t dim_;
	size_t level_;
};

struct AdpativityConfiguration
{
	size_t noRefinements_;
	double threshold_;
	bool type_;
	size_t noPoints_;
	double percent_;
};

struct SLESolverConfiguration
{
	SLESolverType type_;
	double eps_;
	size_t maxIterations_;
	double threshold_;
};

class Learner
{

};

}

}

#endif /* LEARNER_HPP */
