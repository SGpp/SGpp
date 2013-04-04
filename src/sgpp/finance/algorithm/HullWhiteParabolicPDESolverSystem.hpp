/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Chao qi (qic@in.tum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#ifndef HULLWHITEPARABOLICPDESOLVERSYSTEM_HPP
#define HULLWHITEPARABOLICPDESOLVERSYSTEM_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/common/DirichletUpdateVector.hpp"
#include "pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
#include "finance/application/VariableDiscountFactor.hpp"

namespace sg {
  namespace finance {

    /**
     * This class implements the ParabolicPDESolverSystem for the HullWhite
     * Equation.
     */
    class HullWhiteParabolicPDESolverSystem : public sg::pde::OperationParabolicPDESolverSystemFreeBoundaries {
      protected:
        /// theta
        double theta;
        /// sigma
        double sigma;
        /// a
        double a;
        /// the B matrix Operation, on boundary grid
        sg::base::OperationMatrix* OpBBound;
        /// the D matrix Operation, on boundary grid
        sg::base::OperationMatrix* OpDBound;
        /// the E matrix Operation, on boundary grid
        sg::base::OperationMatrix* OpEBound;
        /// the F matrix Operation, on boundary grid
        sg::base::OperationMatrix* OpFBound;
        /// the LTwoDotProduct Operation (Mass Matrix A), on boundary grid
        sg::base::OperationMatrix* OpLTwoBound;

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

        /// vector storing algorithmic dimensions
        std::vector<size_t> HWalgoDims;
        /// Routine to modify the boundaries/inner points of the grid
        sg::base::DirichletUpdateVector* BoundaryUpdate;

        /// access to the variable discount factor
        VariableDiscountFactor* variableDiscountFactor;

        virtual void applyLOperator(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void applyMassMatrix(sg::base::DataVector& alpha, sg::base::DataVector& result);

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param alpha the ansatzfunctions' coefficients
         * @param theta reference to the theta
         * @param sigma reference to the sigma
         * @param a reference to the a
         * @param TimestepSize the size of one timestep used in the ODE Solver
         * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
         *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
         * @param useCoarsen specifies if the grid should be coarsened between timesteps
         * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
         * @param adaptSolveMode adaptivity applied during sloving HullWhite PDE
         * @param numCoarsenPoints number of points that should be coarsened
         * @param refineThreshold surplus threshold for refinement
         * @param refineMode mode used when applying refinements
         * @param refineMaxLevel max. refinement level
         * @param dim_HW dimension of Hull-White (dimension of risk-free rate)
         */
        HullWhiteParabolicPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double sigma, double theta,
                                          double a, double TimestepSize, std::string OperationMode = "ExEul",
                                          bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
                                          int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic",
                                          sg::base::GridIndex::level_type refineMaxLevel = 0, int dim_HW = 1);

        /**
         * Std-Destructor
         */
        virtual ~HullWhiteParabolicPDESolverSystem();

        void finishTimestep();

        /**
         * @param isLastTimestep specify if last time step
         */
        void coarsenAndRefine(bool isLastTimestep = false);

        void startTimestep();

      private:

        /// the dimension of the risk-free rate (Hull-White dimension)
        int dim_r;
    };

  }
}

#endif /* HULLWHITEPARABOLICPDESOLVERSYSTEM_HPP */
