#ifndef LEARNERONLINESGD_HPP
#define LEARNERONLINESGD_HPP

#include <iostream>
#include <cmath>
#include <string>
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/exception/application_exception.hpp"
#include "datadriven/application/Learner.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "datadriven/algorithm/DMSystemMatrix.hpp"


namespace sg {

  namespace datadriven {

    class LearnerOnlineSGD: public sg::datadriven::Learner {

      public:
        LearnerOnlineSGD(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true);

		/*
		 * Implements stochastic gradient descent.
		 *
		 * Note: I cannot pass RefinementFunctor directly because it needs a reference to the alpha DataVector.
		 *
		 * @param trainDataset training dataset: x values
		 * @param classes training dataset: y values
		 * @param GridConfig configuration of initial grid
		 *
		 * @param numIterations number of times SGD is executed before refinement
		 * @param batchSize
		 * @param lambda regularization factor
		 * @param gamma step width
		 *
		 * @param refinement AbstractRefinement
		 * @param refinementType
		 * @param refinementNumPoints
		 *
		 * @param numRuns number of total runs through the dataset (i.e. total number of SGD iterations is numRuns * size of training dataset)
		 * */
		virtual void train(
				sg::base::DataMatrix& trainDataset,
				sg::base::DataVector& classes,
				sg::base::RegularGridConfiguration& GridConfig,

				size_t numIterations,
				size_t batchSize,
				double lambda,
				double gamma,

				sg::base::AbstractRefinement& refinement,
				int refinementCondition,
				int refinementType,
				int refinementNumPoints,

				int numRuns,
				std::ostream **outputStreams_
				);


		virtual ~LearnerOnlineSGD();

		sg::base::DataVector* getAlpha();

      private:
		size_t *batch;
		size_t batchSize;
		size_t batchIndex;
		void pushSGDIndex(size_t index);

		std::ostream **outputStreams;
		void output(int fd, std::string str);

		double getMSE(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes);
		double errorOnMinibatch(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes);
		void performSGDStep(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes, double lambda, double gamma);
		int getRandom(int limit);
	};
  }
}

#endif /* LEARNERSGD_HPP */
