/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)

#ifndef LEARNERSGD_HPP
#define LEARNERSGD_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "datadriven/application/Learner.hpp"

namespace sg {

  namespace datadriven {

    class LearnerSGD: public sg::datadriven::Learner {

      public:
        LearnerSGD(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true);

		/*
		 * Implements stochastic gradient descent.
		 *
		 * @param trainDataset training dataset: x values
		 * @param classes training dataset: y values
		 * @param maxIterations stops after maxIterations
		 * @param eps stop if alpha_i < eps for all i
		 * @param lambda regularization factor
		 * @param gamma step width
		 * */
		virtual void train(
	            sg::base::DataMatrix& trainDataset, 
	            sg::base::DataVector& classes, 
				sg::base::RegularGridConfiguration& GridConfig, 
				size_t maxIterations,
				double eps,
				double lambda,
				double gamma
				);

		virtual ~LearnerSGD();

		sg::base::DataVector* getAlpha();
		sg::base::Grid* getGrid();

      private:
		int getRandom(int limit);
		
	};
  }
}

#endif /* LEARNERSGD_HPP */
