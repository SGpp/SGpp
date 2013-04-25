/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef LEARNERADMMCT_HPP
#define LEARNERADMMCT_HPP

//#include "datadriven/application/LearnerBase.hpp"
#include "parallel/datadriven/application/LearnerAdmm.hpp"
#include "datadriven/tools/TypesDatadriven.hpp"
#include <vector>


using namespace sg::datadriven;

namespace sg {
  namespace parallel {
    class LearnerAdmmCt : public sg::parallel::LearnerAdmm {
      protected:
        /// the grid's coefficients
        sg::base::DataVector* alpha_;
        /// sparse grid object
        sg::base::Grid* grid_;
        // regularization matrix object
        sg::base::OperationMatrix* C_;
        // factorized system matrix
        double* system_matrix;
        /// is verbose output enabled
        bool isVerbose_;
        /// is regression selected
        bool isRegression_;
        /// is the grid trained
        bool isTrained_;
        /// execution time
        double execTime_;
        /// number of executed Floating Point operations
        double GFlop_;
        /// number of transferred Gbytes
        double GByte_;

        int start_index;
        int end_index;
        int grid_size;

        int train_size;

        int admm_iteration;

        int rank;
        int num_tasks;


        double* B_;

        double rho;
        std::vector<int> level_vector;
        double comb_coeff;


        /**
         * abstract method that constructs the corresponding system of linear equations
         * Derived classes MUST overwrite this functions!
         *
         * @param trainDataset training dataset
         * @param lambda lambda regularization parameter
         */
        //sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);

        void createFactorizedSystem(sg::base::DataMatrix& trainDataset, double lambda);

        /**
        * Initialize the grid and its coefficients
        *
        * @param GridCongif structure which describes the regular start grid
        */
        void InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig);

      public:

        /**
         * Constructor
         *
         * @param isRegression
         * @param verbose
         */
        LearnerAdmmCt(const sg::base::RegularGridConfiguration& GridConfig, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true);

        /**
         * Constructor
         *
         * @param tGridFilename path to file that contains a serialized grid
         * @param tAlphaFilenment path to file that contains the grid's coefficients
         * @param isRegression set to true if a regression task should be executed
         * @param verbose set to true in order to allow console output
         */
        //LearnerADMMCT(std::string tGridFilename, std::string tAlphaFilename, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose);

        /**
         * Copy-Constructor
         *
         * @param copyMe LearnerBase that should be duplicated
         */
        //LearnerADMM(const LearnerBase& copyMe);

        /**
         * Destructor
         */
        virtual ~LearnerAdmmCt();

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
        LearnerTiming train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
                            const double lambda, const double convergece_threshold);


      private:
        /*
         * Cholesky factorization subroutine
         */
        int dpotrf(char* uplo, int n, double* a, int lda, int& info);

        /*
        * Linear solver using Cholesky factorization
        */
        int dpotrs(const char* uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info);

    };
  }
}

#endif /* LEARNERADMMCT_HPP */
