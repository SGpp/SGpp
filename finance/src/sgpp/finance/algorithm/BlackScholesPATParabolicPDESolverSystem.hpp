// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESPATPARABOLICPDESOLVERSYSTEM_HPP
#define BLACKSCHOLESPATPARABOLICPDESOLVERSYSTEM_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    /**
     * This class implements the ParabolicPDESolverSystem for the BlackScholes
     * Equation. In this case a principal axis transformation is performed in order
     * to improve the use of spatially adaptive grids and to reduce the
     * calculation effort for higher dimensional cases.
     *
     */
    class BlackScholesPATParabolicPDESolverSystem : public SGPP::pde::OperationParabolicPDESolverSystemFreeBoundaries {
      protected:
        /// the Laplace Operation, on boundary grid
        SGPP::base::OperationMatrix* OpLaplaceBound;
        /// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
        SGPP::base::OperationMatrix* OpLTwoBound;
        /// Eigenvalues of the covariance matrix
        SGPP::base::DataVector* lambda;
        /// Eigenvectors of the covariance matrix
        SGPP::base::DataMatrix* eigenvecs;
        /// Pointer to the mu_hat (transformed drifts and correlation, needed for constraint of American options)
        SGPP::base::DataVector* mu_hat;
        /// use coarsening between timesteps in order to reduce gridsize
        bool useCoarsen;
        /// adaptive mode during solving Black Scholes Equation: coarsen, refine, coarsenNrefine
        std::string adaptSolveMode;
        /// number of points the are coarsened in each coarsening-step !CURRENTLY UNUSED PARAMETER!
        int numCoarsenPoints;
        /// Threshold used to decide if a grid point should be deleted
        float_t coarsenThreshold;
        /// Threshold used to decide if a grid point should be refined
        float_t refineThreshold;
        /// refine mode during solving Black Scholes Equation: classic or maxLevel
        std::string refineMode;
        /// maxLevel max. Level of refinement
        SGPP::base::GridIndex::level_type refineMaxLevel;
        /// the algorithmic dimensions used in this system
        std::vector<size_t> BSalgoDims;
        /// Routine to modify the boundaries/inner points of the grid
        SGPP::base::DirichletUpdateVector* BoundaryUpdate;

        virtual void applyLOperator(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        virtual void applyMassMatrix(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param alpha the ansatzfunctions' coefficients
         * @param lambda eigenvalues of the covariance matrix
         * @param eigenvecs reference to the eigenvectors of the co-variance matrix
         * @param mu_hat reference to transformed drifts and correlation, needed for constraint of American options
         * @param TimestepSize the size of one timestep used in the ODE Solver
         * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
         *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
         * @param useCoarsen specifies if the grid should be coarsened between timesteps
         * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
         * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
         * @param numCoarsenPoints number of point that should be coarsened in one coarsening step !CURRENTLY UNUSED PARAMETER!
         * @param refineThreshold Threshold to decide, if a grid point should be refined
         * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
         * @param refineMaxLevel max. level of refinement during solving
         * @param dStrike the option's strike value
         * @param option_type type of option used here
         */
        BlackScholesPATParabolicPDESolverSystem(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& lambda,
                                                SGPP::base::DataMatrix& eigenvecs, SGPP::base::DataVector& mu_hat, float_t TimestepSize, std::string OperationMode,
                                                float_t dStrike, std::string option_type,
                                                bool useCoarsen = false, float_t coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
                                                int numCoarsenPoints = -1, float_t refineThreshold = 0.0, std::string refineMode = "classic", SGPP::base::GridIndex::level_type refineMaxLevel = 0);

        /**
         * Std-Destructor
         */
        virtual ~BlackScholesPATParabolicPDESolverSystem();

        virtual void finishTimestep();

        /*
         * @ param isLastTimestep specify if last timestep
         */
        virtual void coarsenAndRefine(bool isLastTimestep = false);

        virtual void startTimestep();
    };

  }
}

#endif /* BLACKSCHOLESPATPARABOLICPDESOLVERSYSTEM_HPP */