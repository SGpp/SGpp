/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HEATEQUATIONSOLVER_HPP
#define HEATEQUATIONSOLVER_HPP

#include "pde/application/ParabolicPDESolver.hpp"

#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/common/BoundingBox.hpp"

#include "base/tools/StdNormalDistribution.hpp"

#include "base/application/ScreenOutput.hpp"
#include "base/tools/SGppStopwatch.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

namespace sg {
  namespace pde {

    /**
     * This class provides a simple-to-use solver of the multi dimensional
     * Heat Equation on Sparse Grids.
     *
     * The class's aim is, to hide all complex details of solving the
     * Heat Equation on Sparse Grids!
     *
     * @version $HEAD$
     */
    class HeatEquationSolver : public ParabolicPDESolver {
      protected:
        /// the heat coefficient
        double a;
        /// screen object used in this solver
        sg::base::ScreenOutput* myScreen;

      public:
        /**
         * Std-Constructor of the solver
         */
        HeatEquationSolver();

        /**
         * Std-Destructor of the solver
         */
        virtual ~HeatEquationSolver();

        void constructGrid(sg::base::BoundingBox& myBoundingBox, int level);

        virtual void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        virtual void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        virtual void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul = 0);

        /**
         * This method sets the heat coefficient of the regarded material
         *
         * @param a the heat coefficient
         */
        void setHeatCoefficient(double a);

        /**
         * Inits the grid with a smooth heat distribution based on the
         * normal distribution formula
         *
         * @param alpha reference to the coefficient's vector
         * @param mu the exspected value of the normal distribution
         * @param sigma the sigma of the normal distribution
         * @param factor a factor that is used to stretch the function values
         */
        void initGridWithSmoothHeat(sg::base::DataVector& alpha, double mu, double sigma, double factor);

        /**
         * Inits the screen object
         */
        virtual void initScreen();

        /**
         * Routine to export the RHS of the inner system which has to be
         * solved in order to solve the Poisson equation
         *
         * @param alpha the start solution
         * @param tFilename file into which the rhs is written
         * @param timestepsize the size of the timesteps
         */
        void storeInnerRHS(sg::base::DataVector& alpha, std::string tFilename, double timestepsize);

        /**
         * Routine to export the solution of the inner system which
         * has been calculated by Up/Down scheme
         *
         * @param alpha the start solution
         * @param numTimesteps number timesteps
         * @param timestepsize size of timesteps
         * @param maxCGIterations the maximum of interation in the CG solver
         * @param epsilonCG the epsilon used in the C
         * @param tFilename file into which the rhs is written
         */
        void storeInnerSolution(sg::base::DataVector& alpha, size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, std::string tFilename);
    };

  }
}

#endif /* HEATEQUATIONSOLVER_HPP */
