/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)

#ifndef LEARNERSGD_HPP
#define LEARNERSGD_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/Learner.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace datadriven {

    class LearnerSGD: public SGPP::datadriven::Learner {

      public:
        LearnerSGD(SGPP::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true);

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
	            SGPP::base::DataMatrix& trainDataset, 
	            SGPP::base::DataVector& classes, 
				SGPP::base::RegularGridConfiguration& GridConfig, 
				size_t maxIterations,
				double eps,
				double lambda,
				double gamma
				);

		virtual ~LearnerSGD();

		SGPP::base::DataVector* getAlpha();
		SGPP::base::Grid* getGrid();

      private:
		int getRandom(int limit);
		
	};
  }
}

#endif /* LEARNERSGD_HPP */
