/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LEARNER_HPP
#define LEARNER_HPP

#include "datadriven/application/LearnerBase.hpp"

namespace sg
{

namespace datadriven
{

class Learner : public LearnerBase
{
protected:
	/// regularization operator
	sg::base::OperationMatrix* C_;

	virtual sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix trainDataset, double lambda);

public:
	Learner(sg::base::OperationMatrix& regularization, const bool isRegression, const bool isVerbose = true);

	Learner(const std::string tGridFilename, const std::string tAlphaFilename, sg::base::OperationMatrix& regularization,
			const bool isRegression, const bool isVerbose = true);

	virtual ~Learner();
};

}

}

#endif /* LEARNER_HPP */
