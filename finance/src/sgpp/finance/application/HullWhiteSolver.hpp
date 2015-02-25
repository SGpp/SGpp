// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HULLWHITESOLVER_HPP
#define HULLWHITESOLVER_HPP


#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/application/ScreenOutput.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <sgpp/globaldef.hpp>
#include "../../../../../base/src/sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp"


namespace SGPP {
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


    class HullWhiteSolver : public SGPP::pde::ParabolicPDESolver {
      private:
        ///  the theta value
        float_t theta;
        /// the sigma value
        float_t sigma;
        /// the a value
        float_t a;
        /// the current time
        //float_t t;
        /// stores if the stochastic asset data was passed to the solver
        bool bStochasticDataAlloc;
        /// screen object used in this solver
        SGPP::base::ScreenOutput* myScreen;
        /// use coarsening between timesteps in order to reduce gridsize
        bool useCoarsen;
        /// Threshold used to decide if a grid point should be deleted
        float_t coarsenThreshold;
        /// Threshold used to decide if a grid point should be refined
        float_t refineThreshold;
        /// adaptive mode during solving Black Scholes Equation: none, coarsen, refine, coarsenNrefine
        std::string adaptSolveMode;
        /// refine mode during solving Black Scholes Equation: classic or maxLevel
        std::string refineMode;
        /// number of points the are coarsened in each coarsening-step
        int numCoarsenPoints;
        /// max. level for refinement during solving
        SGPP::base::GridIndex::level_type refineMaxLevel;
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

        void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

        void setStochasticData(float_t theta, float_t sigma, float_t a);


        void solveImplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveExplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, size_t NumImEul = 0);

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
        void initGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType, float_t sigma, float_t a, float_t t, float_t T);

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