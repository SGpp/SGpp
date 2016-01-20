// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESHULLWHITESOLVER_HPP
#define BLACKSCHOLESHULLWHITESOLVER_HPP


#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>


namespace SGPP {
  namespace finance {

    /**
     * This class provides a simple-to-use solver of the multi dimensional Black
     * Scholes Equation that uses Sparse Grids.
     *
     * The class's aim is, to hide all complex details of solving the Black Scholes
     * Equation with Sparse Grids!
     *
     */
    class BlackScholesHullWhiteSolver : public SGPP::pde::ParabolicPDESolver {
      private:
        /// vector that contains the assets' weight
        SGPP::base::DataVector* mus;
        /// vector that contains the standard deviations
        SGPP::base::DataVector* sigmas;
        /// Matrix that contains the correlations
        SGPP::base::DataMatrix* rhos;
        /// the riskfree rate
        float_t r;
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
        /// identifies if the Black Scholes Equation should be solved on a log-transformed grid
        bool useLogTransform;
        /// max. level for refinement during solving
        SGPP::base::GridIndex::level_type refineMaxLevel;
        /// variable to store needed solving iterations
        size_t nNeededIterations;
        /// variable to store the solving time
        float_t dNeededTime;
        /// variable to store start grid size (Inner SGPP::base::Grid)
        size_t staInnerGridSize;
        /// variable to store final grid size (Inner SGPP::base::Grid)
        size_t finInnerGridSize;
        /// variable to store average grid size (Inner SGPP::base::Grid)
        size_t avgInnerGridSize;
        /// Percent how many of the removable points should be tested for deletion
        float_t coarsenPercent;
        /// denotes the number of coarsening procedures within one timestep
        size_t numExecCoarsen;
        float_t theta;
        /// the standard deviations
        float_t sigma;
        ///
        float_t a;
        /// dimension in which BS is calculated
        int dim_BS;
        /// dimension in which HW is calculated
        int dim_HW;

      public:
        /**
         * Std-Constructor of the solver
         */
        BlackScholesHullWhiteSolver(bool useLogTransform = false);

        /**
         * Std-Destructor of the solver
         */
        virtual ~BlackScholesHullWhiteSolver();

        void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

        /**
         * In order to combine the Black Scholes Equation with the Hull White Equation you have to provided
         * some statistical data about the underlying (assets' weight, standard deviation
         * and the correlation between them). This function allows you to set this data for Black-Scholes.
         *
         * @param mus a SGPP::base::DataVector that contains the underlyings' weight
         * @param sigmas a SGPP::base::DataVector that contains the underlyings' standard deviations
         * @param rhos a SGPP::base::DataMatrix that contains the correlations between the underlyings
         * @param r the riskfree rate used in the market model
         * @param theta the theta of HullWhite PDE
         * @param sigma the sigma of HullWhite PDE (vola)
         * @param a the a of HullWhite PDE (mean reversion rate)
         */
        void setStochasticData(SGPP::base::DataVector& mus, SGPP::base::DataVector& sigmas, SGPP::base::DataMatrix& rhos, float_t r, float_t theta, float_t sigma, float_t a);

        /**
         *  defines the dimension of the stoch. processes (BS and HW). default is BS:0, HW:1
         *
         *  @param dim_BS dimension in which Black-Scholes should be calculated
         *  @param dim_HW dimension in which Hull-White should be calculated
         */
        void setProcessDimensions(int dim_BS, int dim_HW);

        void solveImplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveExplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, size_t NumImEul = 0);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         * @param a is the mean reversion rate
         * @param sigma is the volatility
         */
        void initGridWithPayoffBSHW(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType, float_t a, float_t sigma);

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

        /**
         *  enables coarsening of grid during solving the Black Scholes
         *  Equation. The coarsening settings have to be specified in order to
         *  enable coarsening.
         *
         *  @param coarsenThreshold Threshold needed to determine if a grid point should be removed
         *  @param refineMode the Mode used for refining the grid: classic or maxLevel
         *  @param refineMaxLevel max. level for refinement during solving
         *  @param adaptSolveMode adaptive mode during solving equation: coarsen, refine, coarsenNrefine
         *  @param numCoarsenPoints number of points coarsened, -1 all coarsenable points are coarsened
         *  @param refineThreshold Threshold needed to determine if a grid point should be refined
         */
        void setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, float_t coarsenThreshold, float_t refineThreshold);

        /**
         * gets the number of gridpoints at the money
         *
         * Only on Cartesian grids!
         *
         * @param payoffType the payoff type
         * @param strike the option's strike
         * @param eps epsilon to determine the gridpoints, use if at the money is not exactly on grid
         */
        size_t getGridPointsAtMoney(std::string payoffType, float_t strike, float_t eps = 0.0);

        /**
         * gets the number needed iterations to solve Black Scholes Equation
         *
         * @return number of iterations needed to solve Black Scholes Equation, if called before solving 0 is returned
         */
        size_t getNeededIterationsToSolve();

        /**
         * gets needed time in seconds to solve Black Scholes Equation
         *
         * @return needed time in seconds to solve Black Scholes Equation, if called before solving 0 is returned
         */
        float_t getNeededTimeToSolve();

        /**
         * gets the number of points in start grid
         *
         * @returns the number of points in start grid, if called before constructing grid, 0 is returned
         */
        size_t getStartInnerGridSize();

        /**
         * gets the number of points in final grid
         *
         * @returns the number of points in final grid, if called before solving, 0 is returned
         */
        size_t getFinalInnerGridSize();

        /**
         * gets the number of average gridpoints
         *
         * @returns the number of average gridpoints, if called before solving, 0 is returned
         */
        size_t getAverageInnerGridSize();
    };

  }
}

#endif /* BLACKSCHOLESHULLWHITESOLVER_HPP */
