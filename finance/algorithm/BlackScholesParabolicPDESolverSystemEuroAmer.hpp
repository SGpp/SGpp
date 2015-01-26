/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMER_HPP
#define BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMER_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"
#include "finance/tools/Hedging.hpp"

//#define HEDGE
#define HEDGE_EPS 0.05
#define HEDGE_WIDTH_PERCENT 0.95
#define HEDGE_POINTS_PER_DIM 75

namespace sg {
  namespace finance {

    /**
     * This class implements the ParabolicPDESolverSystem for the BlackScholes
     * Equation.
     *
     * Here European or American Options with fix Dirichlet boundaries are solved.
     */
    class BlackScholesParabolicPDESolverSystemEuroAmer : public sg::pde::OperationParabolicPDESolverSystemDirichlet {
      protected:
        /// the riskfree interest rate
        double r;
        /// the delta Operation, on boundary grid
        sg::base::OperationMatrix* OpDeltaBound;
        /// the Gamma Operation, on boundary grid
        sg::base::OperationMatrix* OpGammaBound;
        /// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
        sg::base::OperationMatrix* OpLTwoBound;
        /// the delta Operation, on Inner grid
        sg::base::OperationMatrix* OpDeltaInner;
        /// the Gamma Operation, on Inner grid
        sg::base::OperationMatrix* OpGammaInner;
        /// the LTwoDotProduct Operation (Mass Matrix), on Inner grid
        sg::base::OperationMatrix* OpLTwoInner;
        /// Pointer to the mus
        sg::base::DataVector* mus;
        /// Pointer to the sigmas
        sg::base::DataVector* sigmas;
        /// Pointer to the rhos;
        sg::base::DataMatrix* rhos;
        /// Pointer to the coefficients of operation Delta
        sg::base::DataVector* deltaCoef;
        /// Pointer to the coefficients ot operation Gamma
        sg::base::DataMatrix* gammaCoef;
        /// use coarsening between timesteps in order to reduce gridsize
        bool useCoarsen;
        /// adaptive mode during solving Black Scholes Equation: coarsen, refine, coarsenNrefine
        std::string adaptSolveMode;
        /// number of points the are coarsened in each coarsening-step !CURRENTLY UNUSED PARAMETER!
        int numCoarsenPoints;
        /// Threshold used to decide if a grid point should be deleted
        double coarsenThreshold;
        /// Threshold used to decide if a grid point should be refined
        double refineThreshold;
        /// refine mode during solving Black Scholes Equation: classic or maxLevel
        std::string refineMode;
        /// maxLevel max. Level of refinement
        sg::base::GridIndex::level_type refineMaxLevel;
        /// the algorithmic dimensions used in this system
        std::vector<size_t> BSalgoDims;
        /// store number of executed timesteps
        size_t nExecTimesteps;
        /// the strike of the current option
        double dStrike;
        /// the type of the current option
        std::string option_type;
        /// store whether log coordinates are used
        bool b_log_transform;
#ifdef HEDGE
        /// hedging calculator
        sg::finance::Hedging* myHedge;
#endif

        virtual void applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Build the coefficients for the Gamma Operation, which
         * are the assets' covariance matrix multiplied by 0.5
         *
         * this routine handles also the symmtrie of the
         * gamma operation
         */
        void buildGammaCoefficients();

        /**
         * Build the coefficients for the combined Delta Operation
         */
        void buildDeltaCoefficients();

        /**
         * Build the coefficients for the Gamma Operation, which
         * are the assets' covariance matrix multiplied by 0.5
         *
         * this routine handles also the symmtrie of the
         * gamma operation
         *
         * This function builds the coefficients for the Log Transformed Black Scholes Equation
         */
        void buildGammaCoefficientsLogTransform();

        /**
         * Build the coefficients for the combined Delta Operation
         *
         * This function builds the coefficients for the Log Transformed Black Scholes Equation
         */
        void buildDeltaCoefficientsLogTransform();

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
         * @param dStrike the strike of the current options
         * @param option_type the type of the current option
         * @param bLogTransform indicates that this system belongs to a log-transformed Black Scholes Equation
         * @param useCoarsen specifies if the grid should be coarsened between timesteps
         * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
         * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
         * @param numCoarsenPoints number of point that should be coarsened in one coarsening step !CURRENTLY UNUSED PARAMETER!
         * @param refineThreshold Threshold to decide, if a grid point should be refined
         * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
         * @param refineMaxLevel max. level of refinement during solving
         */
        BlackScholesParabolicPDESolverSystemEuroAmer(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& mu, sg::base::DataVector& sigma,
            sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
            double dStrike, std::string option_type,
            bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
            int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic", sg::base::GridIndex::level_type refineMaxLevel = 0);

        /**
         * Std-Destructor
         */
        virtual ~BlackScholesParabolicPDESolverSystemEuroAmer();

        virtual void finishTimestep();

        virtual void coarsenAndRefine(bool isLastTimestep = false);

        void startTimestep();
    };

  }

}

#endif /* BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMER_HPP */
