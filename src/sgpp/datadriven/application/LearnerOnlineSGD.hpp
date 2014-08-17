/* ****************************************************************************
 * Copyright (C) 2014 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)
#ifndef LEARNERONLINESGD_HPP
#define LEARNERONLINESGD_HPP

#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"
#include "base/exception/application_exception.hpp"
/*
 #include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
 #include "base/operation/BaseOpFactory.hpp"
 #include "base/grid/Grid.hpp"
 #include "base/datatypes/DataVector.hpp"
 #include "base/datatypes/DataMatrix.hpp"
 #include "base/grid/generation/hashmap/AbstractRefinement.hpp"

 #include "datadriven/application/Learner.hpp"
 #include "solver/sle/ConjugateGradients.hpp"
 #include "base/operation/OperationMatrix.hpp"
 #include "datadriven/algorithm/DMSystemMatrix.hpp"
 */

namespace sg {

namespace datadriven {

struct LearnerOnlineSGDRefinementConfiguration {
	std::string refinementCondition;
	size_t numIterations;
	size_t numMinibatchError;

	std::string refinementType;
	int refinementNumPoints;
};

class LearnerOnlineSGD: public sg::datadriven::Learner {

public:
	LearnerOnlineSGD(sg::datadriven::LearnerRegularizationType& regularization,
			const bool isRegression, const bool isVerbose = true);

	virtual void train(sg::base::DataMatrix& mainTrainDataset_,
			sg::base::DataVector& mainClasses_,

			sg::base::DataMatrix& testTrainDataset_,
			sg::base::DataVector& testClasses_,

			sg::base::RegularGridConfiguration& GridConfig,
			sg::datadriven::LearnerOnlineSGDRefinementConfiguration& RefineConfig,
			sg::base::AbstractRefinement& refinement,

			size_t batchSize_, double lambda_, double gamma_,

			int numRuns, std::string errorType_, std::string experimentDir);

	virtual ~LearnerOnlineSGD();

	sg::base::DataVector* getAlpha();

private:
	sg::base::DataMatrix* mainTrainDataset;
	sg::base::DataVector* mainClasses;

	sg::base::DataMatrix* testTrainDataset;
	sg::base::DataVector* testClasses;

	sg::base::DataMatrix* minibatchTrainDataset;
	sg::base::DataVector* minibatchClasses;
	std::list<double>* errorOnMinibatch;

	sg::base::DataVector* errorMB;
	sg::base::DataVector* errorTrainData;

	void pushMinibatch(sg::base::DataVector& x, double y);

	std::string errorType;

	size_t SGDCurrentIndex;
	size_t minibatchSize;
	std::vector<size_t> SGDIndexOrder;
	double lambda;
	double gamma;
	void performSGDStep();

	double getError(sg::base::DataMatrix* trainDataset,
			sg::base::DataVector* classes, sg::base::DataVector *error);
	double getError(std::vector<sg::base::DataVector>* trainDataset,
	    std::vector<double>* classes, sg::base::DataVector *error);
	double getMainError();
	double getTestError();
	double getMinibatchError();
};
}
}

#endif /* LEARNERSGD_HPP */
