/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Chao qi (qic@in.tum.de)

#ifndef HULLWHITESOLVER_HPP
#define HULLWHITESOLVER_HPP


#include "pde/application/ParabolicPDESolver.hpp"

#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/common/BoundingBox.hpp"
#include "base/application/ScreenOutput.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace sg {
  namespace finance {

    /**
     * This class provides a simple-to-use solver of the "multi" dimensional Hull
     * White Equation that uses Sparse Grids.
     *
     * The class's aim is, to hide all complex details of solving the Hull White
     * Equation with Sparse Grids!
     *
     * @version $HEAD$
     */


    class HullWhiteSolver : public sg::pde::ParabolicPDESolver {
      private:
        ///  the theta value
        double theta;
        /// the sigma value
        double sigma;
        /// the a value
        double a;
        /// the current time
        //double t;
        /// stores if the stochastic asset data was passed to the solver
        bool bStochasticDataAlloc;
        /// screen object used in this solver
        sg::base::ScreenOutput* myScreen;
        /// use coarsening between timesteps in order to reduce gridsize
        bool useCoarsen;
        /// Threshold used to decide if a grid point should be deleted
        double coarsenThreshold;
        /// Threshold used to decide if a grid point should be refined
        double refineThreshold;
        /// adaptive mode during solving Black Scholes Equation: none, coarsen, refine, coarsenNrefine
        std::string adaptSolveMode;
        /// refine mode during solving Black Scholes Equation: classic or maxLevel
        std::string refineMode;
        /// number of points the are coarsened in each coarsening-step
        int numCoarsenPoints;
        /// max. level for refinement during solving
        sg::base::GridIndex::level_type refineMaxLevel;
        /// variable to store needed solving iterations


      public:
        /**
         * Std-Constructor of the solver
         */
        HullWhiteSolver();

        /**
         * Std-Destructor of the solver
         */
        virtual ~HullWhiteSolver();

        void constructGrid(sg::base::BoundingBox& myBoundingBox, int level);

        void setStochasticData(double theta, double sigma, double a);


        void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul = 0);

        /**
         * Inits the alpha vector with a payoff function of an European call option
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         * @param sigma the sigma value in HullWhite model
         * @param a the value of a in HullWhite model
         * @param t the current time
         * @param T the maturity time
         */
        void initGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType, double sigma, double a, double t, double T);

        /**
         * Inits the screen object
         */
        void initScreen();

        /**
         * returns the algorithmic dimensions (the dimensions in which the Up Down
         * operations (need for space discretization) should be applied)
         *
         * @return the algorithmic dimensions
         */
        std::vector<size_t> getAlgorithmicDimensions();

        /**
         * sets the algorithmic dimensions (the dimensions in which the Up Down
         * operations (need for space discretization) should be applied)
         *
         * @param newAlgoDims std::vector containing the algorithmic dimensions
         */
        void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

    };

  }
}

#endif /* BLACKSCHOLESSOLVER_HPP */
