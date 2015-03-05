// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODIFIEDBLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP
#define MODIFIEDBLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp>
#include <sgpp/finance/tools/VariableDiscountFactor.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {
    /**
     * This class implements the Modified ParabolicPDESolverSystem for the BlackScholes
     * Equation just use for combination of BlackScholes and HullWhite.
     */
    class ModifiedBlackScholesParabolicPDESolverSystem  : public BlackScholesParabolicPDESolverSystem {
      protected:

        SGPP::base::OperationMatrix* OpFBound;

        virtual void applyLOperator(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

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
         * @param bLogTransform indicates that this system belongs to a log-transformed Black Scholes Equation
         * @param useCoarsen specifies if the grid should be coarsened between timesteps
         * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
         * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
         * @param numCoarsenPoints number of point that should be coarsened in one coarsening step !CURRENTLY UNUSED PARAMETER!
         * @param refineThreshold Threshold to decide, if a grid point should be refined
         * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
         * @param refineMaxLevel max. level of refinement during solving
         * @param dim_HW of Hull-White (= where r value is taken)
         */
        ModifiedBlackScholesParabolicPDESolverSystem(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& mu,
            SGPP::base::DataVector& sigma, SGPP::base::DataMatrix& rho, float_t r, float_t TimestepSize, std::string OperationMode,
            bool bLogTransform, bool useCoarsen, float_t coarsenThreshold, std::string adaptSolveMode,
            int numCoarsenPoints, float_t refineThreshold, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel,
            int dim_HW);

        /**
        * Multiplies the corresponding r coordinates with the whole grid value
        *
        * @param updateVector the vector that should be updated
        */
        virtual void multiplyrBSHW(SGPP::base::DataVector& updateVector);

        /**
        * Std-Destructor
        */
        virtual ~ModifiedBlackScholesParabolicPDESolverSystem();

        virtual void finishTimestep();

        virtual void coarsenAndRefine(bool isLastTimestep = false);

        virtual void startTimestep();

      protected:

        /// the dimension of the risk-free rate (Hull-White dimension)
        int dim_r;

        /// access to the variable discount factor
        VariableDiscountFactor* variableDiscountFactor;
    };

  }
}

#endif /* MODIFIEDBLACKSCHOLESParabolicPDESolverSystem_HPP */