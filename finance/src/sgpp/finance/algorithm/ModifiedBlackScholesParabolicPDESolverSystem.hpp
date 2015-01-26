/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Chao qi (qic@in.tum.de)

#ifndef MODIFIEDBLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP
#define MODIFIEDBLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
#include "finance/tools/VariableDiscountFactor.hpp"

namespace sg {
  namespace finance {
    /**
     * This class implements the Modified ParabolicPDESolverSystem for the BlackScholes
     * Equation just use for combination of BlackScholes and HullWhite.
     */
    class ModifiedBlackScholesParabolicPDESolverSystem  : public BlackScholesParabolicPDESolverSystem {
      protected:

        sg::base::OperationMatrix* OpFBound;

        virtual void applyLOperator(sg::base::DataVector& alpha, sg::base::DataVector& result);

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
        ModifiedBlackScholesParabolicPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& mu,
            sg::base::DataVector& sigma, sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
            bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
            int numCoarsenPoints, double refineThreshold, std::string refineMode, sg::base::GridIndex::level_type refineMaxLevel,
            int dim_HW);

        /**
        * Multiplies the corresponding r coordinates with the whole grid value
        *
        * @param updateVector the vector that should be updated
        */
        virtual void multiplyrBSHW(sg::base::DataVector& updateVector);

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
