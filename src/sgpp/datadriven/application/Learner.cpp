/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/application/Learner.hpp"
#include "datadriven/algorithm/DMSystemMatrix.hpp"

namespace sg
{

namespace datadriven
{

Learner::Learner(sg::base::OperationMatrix& regularization, const bool isRegression, const bool verbose)
	: LearnerBase(isRegression, verbose), C_(&regularization)
{
}

Learner::Learner(const std::string tGridFilename, const std::string tAlphaFilename, sg::base::OperationMatrix& regularization,
		const bool isRegression, const bool verbose)
	: LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose), C_(&regularization)
{
}

Learner::~Learner()
{
}

sg::datadriven::DMSystemMatrixBase* Learner::createDMSystem(sg::base::DataMatrix trainDataset, double lambda)
{
	if (this->grid_ == NULL)
		return NULL;

	return new sg::datadriven::DMSystemMatrix(*(this->grid_), trainDataset, *C_, lambda);
}

}

}
