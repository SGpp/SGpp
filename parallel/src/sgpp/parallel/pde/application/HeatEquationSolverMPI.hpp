// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HEATEQUATIONSOLVERMPI_HPP
#define HEATEQUATIONSOLVERMPI_HPP

#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include <sgpp/globaldef.hpp>
#include "../../../../../../base/src/sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp"


namespace SGPP {
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
    class HeatEquationSolverMPI : public SGPP::pde::ParabolicPDESolver {
      private:
        /// the heat coefficient
        double a;
        /// screen object used in this solver
        SGPP::base::ScreenOutput* myScreen;

      public:
        /**
         * Std-Constructor of the solver
         */
        HeatEquationSolverMPI();

        /**
         * Std-Destructor of the solver
         */
        virtual ~HeatEquationSolverMPI();

        void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

        void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, SGPP::base::DataVector& alpha, size_t NumImEul = 0);

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
        void initGridWithSmoothHeat(SGPP::base::DataVector& alpha, double mu, double sigma, double factor);

        /**
         * Inits the screen object
         */
        void initScreen();
    };

  }
}

#endif /* HEATEQUATIONSOLVERMPI_HPP */