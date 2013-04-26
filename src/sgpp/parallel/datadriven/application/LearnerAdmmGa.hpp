/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef LEARNERBLOCKGAUSSSEIDEL_HPP
#define LEARNERBLOCKGAUSSSEIDEL_HPP

//#include "datadriven/application/LearnerBase.hpp"
#include "parallel/datadriven/application/LearnerADMM.hpp"

#include "datadriven/tools/TypesDatadriven.hpp"

using namespace sg::datadriven;

namespace sg {
  namespace parallel {
    class LearnerADMMGa : public LearnerAdmm {
      public:
        /**
         * Constructor
         *
         * @param isRegression
         * @param verbose
         */
        LearnerADMMGa(const sg::base::RegularGridConfiguration& GridConfig, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true, SolutionType method = Cholesky)
          : LearnerAdmm(GridConfig, regularization, isRegression, isVerbose, method) {};

        LearnerADMMGa(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true, SolutionType method = Cholesky)
          : LearnerAdmm(regularization, isRegression, isVerbose, method) {};
        /**
         * Constructor
         *
         * @param tGridFilename path to file that contains a serialized grid
         * @param tAlphaFilenment path to file that contains the grid's coefficients
         * @param isRegression set to true if a regression task should be executed
         * @param verbose set to true in order to allow console output
         */
        LearnerADMMGa(std::string tGridFilename, std::string tAlphaFilename, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose, SolutionType method = Cholesky)
          : LearnerAdmm(tGridFilename, tAlphaFilename, regularization, isRegression, isVerbose, method) {};

        /**
         * Learning a dataset with spatially adaptive sparse grids
         *
         * @param testDataset the training dataset
         * @param classes classes corresponding to the training dataset
         * @param GridConfig configuration of the regular start grid
         * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
         * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
         * @param AdaptConfig configuration of the adaptivity strategy
         * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
         * @param lambda regularization parameter lambda
         */
        virtual LearnerTiming train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
                                    const double lambda, const double convergece_threshold, double processes = -1);

    };
  }
}

#endif /* LEARNERBLOCKGAUSSSEIDEL_HPP */
