/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author author Chao qi (qic@in.tum.de)

#ifndef BLACKSCHOLESHULLWHITESOLVER_HPP
#define BLACKSCHOLESHULLWHITESOLVER_HPP


#include "pde/application/ParabolicPDESolver.hpp"

#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/common/BoundingBox.hpp"

#include "base/tools/StdNormalDistribution.hpp"

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
     * This class provides a simple-to-use solver of the multi dimensional Black
     * Scholes Equation that uses Sparse Grids.
     *
     * The class's aim is, to hide all complex details of solving the Black Scholes
     * Equation with Sparse Grids!
     *
     * @version $HEAD$
     */
    class BlackScholesHullWhiteSolver : public sg::pde::ParabolicPDESolver {
      private:
        /// vector that contains the assets' weight
        sg::base::DataVector* mus;
        /// vector that contains the standard deviations
        sg::base::DataVector* sigmas;
        /// Matrix that contains the correlations
        sg::base::DataMatrix* rhos;
        /// the riskfree rate
        double r;
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
        /// identifies if the Black Scholes Equation should be solved on a log-transformed grid
        bool useLogTransform;
        /// max. level for refinement during solving
        sg::base::GridIndex::level_type refineMaxLevel;
        /// variable to store needed solving iterations
        size_t nNeededIterations;
        /// variable to store the solving time
        double dNeededTime;
        /// variable to store start grid size (Inner sg::base::Grid)
        size_t staInnerGridSize;
        /// variable to store final grid size (Inner sg::base::Grid)
        size_t finInnerGridSize;
        /// variable to store average grid size (Inner sg::base::Grid)
        size_t avgInnerGridSize;
        /// Percent how many of the removable points should be tested for deletion
        double coarsenPercent;
        /// denotes the number of coarsening procedures within one timestep
        size_t numExecCoarsen;
        double theta;
        /// the standard deviations
        double sigma;
        ///
        double a;
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

        void constructGrid(sg::base::BoundingBox& myBoundingBox, int level);

        /**
         * In order to combine the Black Scholes Equation with the Hull White Equation you have to provided
         * some statistical data about the underlying (assets' weight, standard deviation
         * and the correlation between them). This function allows you to set this data for Black-Scholes.
         *
         * @param mus a sg::base::DataVector that contains the underlyings' weight
         * @param sigmas a sg::base::DataVector that contains the underlyings' standard deviations
         * @param rhos a sg::base::DataMatrix that contains the correlations between the underlyings
         * @param r the riskfree rate used in the market model
         * @param theta the theta of HullWhite PDE
         * @param sigma the sigma of HullWhite PDE (vola)
         * @param a the a of HullWhite PDE (mean reversion rate)
         */
        void setStochasticData(sg::base::DataVector& mus, sg::base::DataVector& sigmas, sg::base::DataMatrix& rhos, double r, double theta, double sigma, double a);

        /**
         *  defines the dimension of the stoch. processes (BS and HW). default is BS:0, HW:1
         *
         *  @param dim_BS dimension in which Black-Scholes should be calculated
         *  @param dim_HW dimension in which Hull-White should be calculated
         */
        void setProcessDimensions(int dim_BS, int dim_HW);

        void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul = 0);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         * @param a is the mean reversion rate
         * @param sigma is the volatility
         */
        void initGridWithPayoffBSHW(sg::base::DataVector& alpha, double strike, std::string payoffType, double a, double sigma);

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
        void setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, sg::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, double coarsenThreshold, double refineThreshold);

        /**
         * gets the number of gridpoints at the money
         *
         * Only on Cartesian grids!
         *
         * @param payoffType the payoff type
         * @param strike the option's strike
         * @param eps epsilon to determine the gridpoints, use if at the money is not exactly on grid
         */
        size_t getGridPointsAtMoney(std::string payoffType, double strike, double eps = 0.0);

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
        double getNeededTimeToSolve();

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
