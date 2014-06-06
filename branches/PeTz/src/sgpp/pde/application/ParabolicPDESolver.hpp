/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef PARABOLICPDESOLVER_HPP
#define PARABOLICPDESOLVER_HPP

#include "pde/application/PDESolver.hpp"

namespace sg {
  namespace pde {

    /**
     * This class extends the PDESolver with functions that are needed to
     * solve parabolic PDEs
     *
     * @version $HEAD$
     */
    class ParabolicPDESolver : public PDESolver {
      protected:
        /// the size of one timestep
        //double timestepSize;
        /// The number of timesteps that are executed during solving
        //size_t nTimesteps;

      public:
        /**
         * Std-Constructor of the solver
         */
        ParabolicPDESolver();

        /**
         * Std-Destructor of the solver
         */
        virtual ~ParabolicPDESolver();

        /**
         * Call this routine to use an explicit Euler algorithm to solve the parabolic PDE
         *
         * @param numTimesteps the number of timesteps that should be executed
         * @param timestepsize the size of the interval one timestep moves forward
         * @param maxCGIterations the maximum of interation in the CG solver
         * @param epsilonCG the epsilon used in the CG
         * @param alpha the coefficients of the Sparse Gird's basis functions
         * @param verbose enables verbose output during solving
         * @param generateAnimation set this to true, if you want to generate a grid output in every timestep
         * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
         */
        virtual void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20) = 0;

        /**
         * Call this routine to use an explicit Euler algorithm to solve the parabolic PDE
         *
         * @param numTimesteps the number of timesteps that should be executed
         * @param timestepsize the size of the interval one timestep moves forward
         * @param maxCGIterations the maximum of interation in the CG solver
         * @param epsilonCG the epsilon used in the CG
         * @param alpha the coefficients of the Sparse Gird's basis functions
         * @param verbose enables verbose output during solving
         * @param generateAnimation set this to true, if you want to generate a grid output in every timestep
         * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
         */
        virtual void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20) = 0;

        /**
         * Call this routine to use the Crank Nicolson algorithm to solve the parabolic PDE
         *
         * @param numTimesteps the number of timesteps that should be executed
         * @param timestepsize the size of the interval one timestep moves forward
         * @param maxCGIterations the maximum of interation in the CG solver
         * @param epsilonCG the epsilon used in the CG
         * @param alpha the coefficients of the Sparse Gird's basis functions
         * @param NumImEul specifies how many ImEul steps should be executed before CrNic is used, default is 0
         */
        virtual void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul = 0) = 0;
    };

  }
}

#endif /* PARABOLICPDESOLVER_HPP */
