/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HEATEQUATIONPARABOLICPDESOLVERSYSTEM_HPP
#define HEATEQUATIONPARABOLICPDESOLVERSYSTEM_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"
#include "pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"

namespace sg {
  namespace pde {

    /**
     * This class implements the ParabolicPDESolverSystem for the
     * Heat Equation.
     */
    class HeatEquationParabolicPDESolverSystem : public OperationParabolicPDESolverSystemDirichlet {
      private:
        /// the heat coefficient
        double a;
        /// the Laplace Operation (Stiffness Matrix), on boundary grid
        sg::base::OperationMatrix* OpLaplaceBound;
        /// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
        sg::base::OperationMatrix* OpMassBound;
        /// the Laplace Operation (Stiffness Matrix), on inner grid
        sg::base::OperationMatrix* OpLaplaceInner;
        /// the LTwoDotProduct Operation (Mass Matrix), on inner grid
        sg::base::OperationMatrix* OpMassInner;

        void applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

        void applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

        void applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result);

        void applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result);

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param alpha the sparse grid's coefficients
         * @param a the heat coefficient
         * @param TimestepSize the size of one timestep used in the ODE Solver
         * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
         *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
         */
        HeatEquationParabolicPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double a, double TimestepSize, std::string OperationMode = "ExEul");

        /**
         * Std-Destructor
         */
        virtual ~HeatEquationParabolicPDESolverSystem();

        virtual void finishTimestep();

        virtual void coarsenAndRefine(bool isLastTimestep = false);

        virtual void startTimestep();
    };

  }
}

#endif /* HEATEQUATIONPARABOLICPDESOLVERSYSTEM_HPP */
