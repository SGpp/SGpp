/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELOMP_HPP
#define BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELOMP_HPP

#include "finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"

namespace sg {
  namespace finance {

    /**
     * This class implements the ParabolicPDESolverSystem for the BlackScholes
     * Equation.
     *
     *
     * Here European or American Options with fix Dirichlet boundaries are solved.
     *
     * It's derived from the existing BlackScholesParabolicPDESolverSystemEuropean but uses
     * the OMP task concept to enable further parallelization possibilities
     * in the calculation of the space-discretization operator (L)
     */
    class BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP : public BlackScholesParabolicPDESolverSystemEuroAmer {
      protected:
        virtual void applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param alpha the ansatzfunctions' coefficients
         * @param mu reference to the mus
         * @param sigma reference to the sigmas
         * @param rho reference to the rhos
         * @param r the riskfree interest rate
         * @param TimestepSize the size of one timestep used in the ODE Solver
         * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
         *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
         * @param dStrike strike
         * @param option_type type of option
         * @param bLogTransform indicates that this system belongs to a log-transformed Black Scholes Equation
         * @param useCoarsen specifies if the grid should be coarsened between timesteps
         * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
         * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
         * @param numCoarsenPoints number of point that should be coarsened in one coarsening step !CURRENTLY UNUSED PARAMETER!
         * @param refineThreshold Threshold to decide, if a grid point should be refined
         * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
         * @param refineMaxLevel max. level for refinement during solving
         */
        BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& mu, sg::base::DataVector& sigma,
            sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
            double dStrike, std::string option_type,
            bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
            int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic", sg::base::GridIndex::level_type refineMaxLevel = 0);

        /**
         * Std-Destructor
         */
        virtual ~BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP();

        /**
         * Multiplicates a vector with the matrix, parallel
         *
         * @param alpha sg::base::DataVector that contains the ansatzfunctions' coefficients
         * @param result sg::base::DataVector into which the result of the space discretization operation is stored
         */
        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * generates the right hand side of the system, parallel
         *
         * @return returns the rhs
         */
        virtual sg::base::DataVector* generateRHS();
    };

  }

}

#endif /* BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELOMP_HPP */
