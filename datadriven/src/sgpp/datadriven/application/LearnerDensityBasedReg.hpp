// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERDENSITYBASEDREG_HPP_
#define LEARNERDENSITYBASEDREG_HPP_

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/pde/application/RegularizationConfiguration.hpp>

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Class that implements a regression learner
     * using density estimation on Sparse Grids.
     *
     */
    class LearnerDensityBasedReg: public LearnerBase {
      public:

        /**
         * Constructor
         *
         * @param regularization regularization type
         * @param border offset for the normalization of the class data
         */
        LearnerDensityBasedReg(
          SGPP::pde::RegularizationType& regularization,
          float_t border = 0.);

        /**
         * Destructor
         */
        virtual ~LearnerDensityBasedReg();

        /**
         * Learning a dataset with spatially adaptive sparse grids
         *
         * @param trainDataset the training dataset
         * @param classes classes corresponding to the training dataset
         * @param GridConfig configuration of the regular start grid
         * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
         * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
         * @param AdaptConfig configuration of the adaptivity strategy
         * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
         * @param lambda regularization parameter lambda
         */
        virtual LearnerTiming train(SGPP::base::DataMatrix& trainDataset,
                                    SGPP::base::DataVector& classes,
                                    const SGPP::base::RegularGridConfiguration& GridConfig,
                                    const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
                                    const SGPP::solver::SLESolverConfiguration& SolverConfigFinal,
                                    const SGPP::base::AdpativityConfiguration& AdaptConfig,
                                    bool testAccDuringAdapt, const float_t lambda);

        /**
         * Executes a regression test for a given dataset and returns the result
         *
         * @param testDataset dataset that is evaluated with the current learner
         *
         * @return regression values of testDataset
         */
        virtual SGPP::base::DataVector predict(SGPP::base::DataMatrix& testDataset);

        /**
         * simple dump of the one dimensional density function for a specific data point
         *  into file, e.g. used to plot with gnuplot.
         *
         * @param point point for which the one dimensional density function is computed
         * @param fileName filename to store the dump to
         * @param resolution resolution of function plot
         */
        void dumpDensityAtPoint(SGPP::base::DataVector& point, std::string fileName,
                                unsigned int resolution);

      protected:
        /// regularization mode
        SGPP::pde::RegularizationType CMode_;
        /// regularization operator
        SGPP::base::OperationMatrix* C_;
        /// maximum value (used for de-normalization)
        float_t maxValue_;
        /// minimum value (used for de-normalization)
        float_t minValue_;
        /// border for normalization of the class vector
        float_t border_;

        /**
         * inherited from LearnerBase, but not used
         */
        virtual SGPP::datadriven::DMSystemMatrixBase* createDMSystem(
          SGPP::base::DataMatrix& trainDataset, float_t lambda);

    };

  } /* namespace datadriven */
} /* namespace SGPP */
#endif /* LEARNERDENSITYBASEDREG_HPP_ */
