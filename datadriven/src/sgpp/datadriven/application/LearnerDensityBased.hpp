// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LEARNERDENSITYBASED_HPP_
#define LEARNERDENSITYBASED_HPP_

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>

#include <sgpp/globaldef.hpp>
#include "../../../../../base/src/sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp"


namespace SGPP {

  namespace datadriven {


    class LearnerDensityBased: public SGPP::datadriven::LearnerBase {
      protected:
        //Mapping from class index to class number:
        std::map<int, double> index_to_class_;
        //Stores the coefficients for every class
        std::vector<SGPP::base::DataVector> alphaVec_;
        /// regularization mode
        SGPP::datadriven::LearnerRegularizationType CMode_;
        //with prior
        bool withPrior;
        //number of classes
        size_t nrClasses;
        // prior of data
        std::vector<double> prior;
        // vectors of grids
        std::vector<SGPP::base::Grid*> gridVec_;
        // vector of regterms
        std::vector<SGPP::base::OperationMatrix*> CVec_;
      public:
        LearnerDensityBased(SGPP::datadriven::LearnerRegularizationType&, const bool isRegression, const bool isVerbose = true);
        virtual ~LearnerDensityBased();

        /**
         * Create a grid for each class
         *
         * @param GridConfig grid config
         */
        virtual void InitializeGrid(const SGPP::base::RegularGridConfiguration& GridConfig);

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
        virtual LearnerTiming train(SGPP::base::DataMatrix& testDataset, SGPP::base::DataVector& classes,
                                    const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
                                    const SGPP::solver::SLESolverConfiguration& SolverConfigFinal, const SGPP::base::AdpativityConfiguration& AdaptConfig,
                                    bool testAccDuringAdapt, const double lambda);


        virtual SGPP::base::DataVector predict(SGPP::base::DataMatrix& testDataset);
        /// construct system matrix
        virtual SGPP::datadriven::DMSystemMatrixBase* createDMSystem(SGPP::base::DataMatrix& trainDataset, double lambda);

        /**
         * Returns the execution time
         */
        time_t getExecTime();

        /**
         * Returns number of grid points for the density
         * with the maximum number of grid points
         */
        size_t getNrGridPoints();

        /**
         * Get Prior
         */
        bool getWithPrior() {
          return withPrior;
        }

        /**
         * Set prior
         *
         * @param p prior
         */
        bool setWithPrior(bool p) {
          withPrior = p;
          return withPrior;
        }

        /**
         * Get number of classes
         */
        size_t getNrClasses() {
          return nrClasses;
        }

        /**
         * Set number of classes
         *
         * @param c set number of classes
         */
        size_t setNrClasses(size_t c) {
          nrClasses = c;
          return nrClasses;
        }



    };

  }

}

#endif /* LEARNERDENSITYBASED_HPP_ */