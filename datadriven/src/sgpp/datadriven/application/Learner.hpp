/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LEARNER_HPP
#define LEARNER_HPP

#include <sgpp/datadriven/application/LearnerBase.hpp>

namespace sg {

  namespace datadriven {

    enum LearnerRegularizationType {
      Laplace,
      Identity
    };

    /**
     * This class implements standard sparse grid regression
     * with an arbitrary regularization operator
     */
    class Learner : public LearnerBase {
      protected:
        /// regularization mode
        sg::datadriven::LearnerRegularizationType CMode_;
        /// regularization operator
        sg::base::OperationMatrix* C_;

        /// construct system matrix
        virtual sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);

      public:
        /**
         * Constructor
         *
         * @param regularization enum that gives the regurlarization method
         * @param isRegression flag for regression
         * @param isVerbose flag for verbose output
         */
        Learner(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true);

        /**
         * Constructor
         *
         * @param tGridFilename path to file that contains a serialized grid
         * @param tAlphaFilename path to file that contains the grid's coefficients
         * @param regularization enum that gives the regurlarization method
         * @param isRegression set to true if a regression task should be executed
         * @param isVerbose set to true in order to allow console output
         */
        Learner(const std::string tGridFilename, const std::string tAlphaFilename, sg::datadriven::LearnerRegularizationType& regularization,
                const bool isRegression, const bool isVerbose = true);

        /**
         * Destructor
         */
        virtual ~Learner();
    };

  }

}

#endif /* LEARNER_HPP */
