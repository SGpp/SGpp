/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Peter Hoffmann (peter.hoffmann@mytum.de)

#ifndef OPERATIONPARABOLICPDESOLVERSYSTEM_HPP
#define OPERATIONPARABOLICPDESOLVERSYSTEM_HPP

#include "base/grid/Grid.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {

    /**
     * Abstract definition of a System that is used to solve parabolic partial
     * differential equations. So an instance of this class has to pass to
     * any ODE Solver used in SGpp.
     *
     * \f$A \dot{u} = L \vec{u}\f$
     *
     * A: mass matrix
     * L: space discretization (L-Operator)
     */
    class OperationParabolicPDESolverSystem : public sg::base::OperationMatrix {
      protected:
        /// Pointer to the alphas (ansatzfunctions' coefficients)
        sg::base::DataVector* alpha_complete;
        /// Pointer to the alphas from the last timestep, needed when using variable timestep sizes
        sg::base::DataVector* alpha_complete_old;
        /// Pointer to temporary alphas, needed when using variable timestep sizes
        sg::base::DataVector* alpha_complete_tmp;
        /// Pointer to the grid from the last iteration
        sg::base::GridStorage* oldGridStorage;
        /// Pointer to the grid from the last aborted iteration
        sg::base::GridStorage* secondGridStorage;

        /**
         *  specifies in which solver this matrix is used, valid values are:
         *  ExEul for explicit Euler
         *  ImEul for implicit Euler
         *  CrNic for Crank Nicolson solver
         */
        std::string tOperationMode;
        /// the size of one timestep used in the ODE Solver
        double TimestepSize;
        /// the size of the last timestep
        double TimestepSize_old;

        sg::base::DataVector* rhs;
        /// Pointer to the grid object
        sg::base::Grid* BoundGrid;
        /// Stores number of average gridpoints, inner grid
        size_t numSumGridpointsInner;
        /// Stores number of average gridpoints, complete grid
        size_t numSumGridpointsComplete;

        /// checks whether a new ODE solver has been selected after creation
        bool bnewODESolver;

      public:
        /**
         * Constructor
         */
        OperationParabolicPDESolverSystem();

        /**
         * Destructor
         */
        virtual ~OperationParabolicPDESolverSystem();

        /**
         * Multiplicates a vector with the matrix
         *
         * @param alpha sg::base::DataVector that contains the ansatzfunctions' coefficients
         * @param result sg::base::DataVector into which the result of the space discretization operation is stored
         */
        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

        /**
         * generates the right hand side of the system
         *
         * @return returns the rhs
         */
        virtual sg::base::DataVector* generateRHS() = 0;

        /**
         * performs some action that might be needed after a timestep has be finished in the ODE
         * Solver, e.g. some boundary adjustments.
         */
        virtual void finishTimestep() = 0;

        virtual void coarsenAndRefine(bool isLastTimestep = false) = 0;


        /**
         * Implements some start jobs of every timestep, e.g.discounting boundaries
         */
        virtual void startTimestep() = 0;

        /**
         * get the pointer to the underlying grid object
         *
         * @return returns a pointer to the underlying grid object
         */
        sg::base::Grid* getGrid();

        /**
         * gets a pointer to the sparse grids coefficients used in the CG method to solve
         * one timestep. This is useful because (direchlet) boundaries can be skipped when
         * solving the system.
         *
         * @return alpha vector for CG method
         */
        virtual sg::base::DataVector* getGridCoefficientsForCG() = 0;

        /**
         * gets a pointer to the sparse grids coefficients with evtl. boundaries
         *
         * @return alpha vector of complete grid
         */
        sg::base::DataVector* getGridCoefficients();

        /**
         * defines the used ODE Solver for this instance, this is important because
         * the implementation of mult and generateRHS depends on the used
         * ODE solver
         *
         * @param ode the used ODESolver: ExEul, ImEul or CrNic
         */
        void setODESolver(std::string ode);

        /**
         * Returns the specified ODE solver for this instance
         *
         * @return the ODE solver: ExEul, ImEul or CrNic
         */
        std::string getODESolver();

        /**
         * Returns the number of average grid points for the complete grid
         *
         * @return the number of average grid points for the complete grid
         */
        size_t getSumGridPointsComplete();

        /**
         * Returns the number of average grid points for the inner grid
         *
         * @return the number of average grid points for the inner grid
         */
        size_t getSumGridPointsInner();

        /**
         * set the size of the new timestep
         *
         * @param newTimestepSize the size of the next timestep
         */
        void setTimestepSize(double newTimestepSize);

        /**
         * aborts the current timestep execution
         */
        void abortTimestep();

        /**
         * stores the current alpha_complete into alpha_complete_old to be available in the next timestep
         */
        void saveAlpha();

        /**
         * stores the values of the (dehierarchized) grid in the sg::base::DataVector Values used by time step size control methods
         *
         * @param Values sg::base::DataVector in which the values will be stored
         */
        void getGridCoefficientsForSC(sg::base::DataVector& Values);

        sg::base::GridStorage* getGridStorage();

        sg::base::GridStorage* getOldGridStorage();

        sg::base::GridStorage* getSecondGridStorage();
    };

  }
}

#endif /* OPERATIONPARABOLICPDESOLVERSYSTEM_HPP */
