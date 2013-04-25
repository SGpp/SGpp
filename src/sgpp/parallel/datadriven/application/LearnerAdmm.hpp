/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef LEARNERADMM_HPP
#define LEARNERADMM_HPP

//#include "datadriven/application/LearnerBase.hpp"
#include "datadriven/application/Learner.hpp"

#include "datadriven/tools/TypesDatadriven.hpp"

using namespace sg::datadriven;

namespace sg {
  namespace parallel {
    class LearnerAdmm : public sg::datadriven::Learner {
      public:

        enum SolutionType {
          Cholesky, CG
        };

        /**
         * Constructor
         *
         * @param isRegression
         * @param verbose
         */
        LearnerAdmm(const sg::base::RegularGridConfiguration& GridConfig, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true, SolutionType method = Cholesky);

        LearnerAdmm(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose = true, SolutionType method = Cholesky);

        /**
         * Constructor
         *
         * @param tGridFilename path to file that contains a serialized grid
         * @param tAlphaFilenment path to file that contains the grid's coefficients
         * @param isRegression set to true if a regression task should be executed
         * @param verbose set to true in order to allow console output
         */
        LearnerAdmm(std::string tGridFilename, std::string tAlphaFilename, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose, SolutionType method = Cholesky);

        /**
         * Copy-Constructor
         *
         * @param copyMe LearnerBase that should be duplicated
         */
        //LearnerADMM(const LearnerBase& copyMe);

        /**
         * Destructor
         */
        virtual ~LearnerAdmm();

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
                                    const double lambda, const double convergece_threshold);

        void set_rho(double rho);
        double get_rho();

        //virtual sg::base::DataVector predict(sg::base::DataMatrix& testDataset);
      protected:
        // factorized system matrix
        double* system_matrix;

        int first_index;
        int last_index;
        int subproblem_length;

        int train_size;

        int admm_iteration;

        int rank;
        int num_tasks;

        double* B_;

        double rho;

        SolutionType solution_method;

        double cg_residual_norm_threshold;


        /**
         * abstract method that constructs the corresponding system of linear equations
         * Derived classes MUST overwrite this functions!
         *
         * @param trainDataset training dataset
         * @param lambda lambda regularization parameter
         */
        //sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);

        virtual void createFactorizedSystem(sg::base::DataMatrix& trainDataset, double lambda);

        virtual void init();

        /**
        * Initialize the grid and its coefficients
        *
        * @param GridCongif structure which describes the regular start grid
        */
        void InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig);

        /**
         * Solve system using CG
         */
        void solve_cg(double* b, double* x, int& info);


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

#endif /* LEARNERADMM_HPP */
