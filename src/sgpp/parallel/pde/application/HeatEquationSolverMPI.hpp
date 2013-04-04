/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HEATEQUATIONSOLVERMPI_HPP
#define HEATEQUATIONSOLVERMPI_HPP

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
  namespace parallel {

    /**
     * This class provides a simple-to-use solver of the multi dimensional
     * Heat Equation on Sparse Grids.
     *
     * The class's aim is, to hide all complex details of solving the
     * Heat Equation on Sparse Grids!
     *
     * This version offers support for MPI parallelization!
     *
     * @version $HEAD$
     */
    class HeatEquationSolverMPI : public sg::pde::ParabolicPDESolver {
      private:
        /// the heat coefficient
        double a;
        /// screen object used in this solver
        sg::base::ScreenOutput* myScreen;

      public:
        /**
         * Std-Constructor of the solver
         */
        HeatEquationSolverMPI();

        /**
         * Std-Destructor of the solver
         */
        virtual ~HeatEquationSolverMPI();

        void constructGrid(sg::base::BoundingBox& myBoundingBox, int level);

        void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul = 0);

        /**
         * This method sets the heat coefficient of the regarded material
         *
         * @param a the heat coefficient
         */
        void setHeatCoefficient(double a);

        /**
         * Inits the grid with an smooth heat distribution based the
         * normal distribution formula
         *
         * @param alpha reference to the coefficients vector
         * @param mu the exspected value of the normal distribution
         * @param sigma the sigma of the normal distribution
         * @param factor a factor that is used to stretch the function values
         */
        void initGridWithSmoothHeat(sg::base::DataVector& alpha, double mu, double sigma, double factor);

        /**
         * Inits the screen object
         */
        void initScreen();
    };

  }
}

#endif /* HEATEQUATIONSOLVERMPI_HPP */
