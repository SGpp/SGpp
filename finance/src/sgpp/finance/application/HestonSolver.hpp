// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef HESTONSOLVER_HPP
#define HESTONSOLVER_HPP


#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/solver/ODESolver.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

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
     * This class provides abstracts the functionality for solving the multi dimensional Heston
     * Equation on Sparse Grids.
     *
     * @version $HEAD$
     */
    class HestonSolver : public SGPP::pde::ParabolicPDESolver {
      protected:

        /// Vector that contains the thetas
        SGPP::base::DataVector* thetas;

        /// Vector that contains the kappas
        SGPP::base::DataVector* kappas;

        /// Vector that contains the volatilities of the volatilities
        SGPP::base::DataVector* volvols;

        /// Matrix that contains the correlations
        SGPP::base::DataMatrix* hMatrix;

        // The number of assets (half the dimension of the PDE)
        size_t numAssets;

        /// The market riskfree rate
        float_t r;

        /// Stores if the stochastic asset data has been passed to the solver
        bool bStochasticDataAlloc;

        /// Screen object used in this solver
        SGPP::base::ScreenOutput* myScreen;

        /// Use coarsening between timesteps in order to reduce gridsize
        bool useCoarsen;

        /// Threshold used to decide if a grid point should be deleted
        float_t coarsenThreshold;

        /// Threshold used to decide if a grid point should be refined
        float_t refineThreshold;

        /// Adaptive mode during solving Heston Equation: none, coarsen, refine, coarsenNrefine
        std::string adaptSolveMode;

        /// Refine mode during solving Heston Equation: classic or maxLevel
        std::string refineMode;

        /// Number of points the are coarsened in each coarsening-step
        int numCoarsenPoints;

        /// Identifies if the Heston Equation should be solved using log-transformed stock-price coordinates. The variance coordinates are always linear.
        bool useLogTransform;

        /// Max. level for refinement during solving
        SGPP::base::GridIndex::level_type refineMaxLevel;

        /// Variable to store needed solving iterations
        size_t nNeededIterations;

        /// Variable to store the solving time
        float_t dNeededTime;

        /// Variable to store start grid size (Inner SGPP::base::Grid)
        size_t staInnerGridSize;

        /// Variable to store final grid size (Inner SGPP::base::Grid)
        size_t finInnerGridSize;

        /// Variable to store average grid size (Inner SGPP::base::Grid)
        size_t avgInnerGridSize;

        /// Type of the Option to solve
        std::string tBoundaryType;

        /// Stores the current time until which the option has been solved
        float_t current_time;

        /// Stores the strike of the current option
        float_t dStrike;

        /// Stores the option type of the current option
        std::string payoffType;

        /**
         * Returns the option value (payoff value) for an European call option
         *
         * @param assetValue the current asset's value
         * @param strike the strike price of the option
         *
         * @return the call premium
         */
        virtual float_t get1DEuroCallPayoffValue(float_t assetValue, float_t strike);

        /**
         * Initialises the alpha vector with a payoff function of an European call option or put option.
         * The grid is initialized based on Cartesian coordinates!
         *
         * @param alpha the coefficient vector of the grid's basis functions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initCartesianGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType);

        /**
         * Initialises the alpha vector with a payoff function of an European call option or put option
         * The grid is initialized based on log-transformed coordinates!
         *
         * @param alpha the coefficient vector of the grid's basis functions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initLogTransformedGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType);

      public:
        /**
         * Std-Constructor of the solver
         *
         * @param useLogTransform speciefies if a log transformed formulation should be used for solving the Heston Equation
         */
        HestonSolver(bool useLogTransform = true);

        /**
         * Std-Destructor of the solver
         */
        virtual ~HestonSolver();

        /**
         * Builds the sparse grid for use in the solver.
         *
         * @param myBoundingBox bounding box for the sparse grid
         * @param level sparse grid level
         */
        virtual void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

        /**
         * This function tries to refine the grid such that
         * most of the grid points are used for interpolation of the singularity. So this grid
         * is able to approximate the start solution better.
         *
         * After refining the grid the payoff function is applied to the grid.
         *
         * Only on Cartesian grids!
         *
         * @param alpha reference to a SGPP::base::DataVector object that contains the gird ansatzfunction's coefficients
         * @param strike containing the option's strike
         * @param payoffType the type of payoff Function used ONLY supported: avgM
         * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
         */
        virtual void refineInitialGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType, float_t dStrikeDistance);

        /**
         * This function tries to refine the grid such that
         * most of the grid points are used for interpolation of the singularity. So this grid
         * is able to approximate the start solution better. Refining is done only if the max
         * refinement level hasn't be reached.
         *
         * After refining the grid the payoff function is applied to the grid.
         *
         * Only on Cartesian grids!
         *
         * @param alpha reference to a SGPP::base::DataVector object that contains the gird ansatzfunction's coefficients
         * @param strike containing the option's strike
         * @param payoffType the type of payoff Function used ONLY supported: avgM
         * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
         * @param maxLevel maximum level of refinement
         */
        virtual void refineInitialGridWithPayoffToMaxLevel(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType, float_t dStrikeDistance, SGPP::base::GridIndex::level_type maxLevel);

        /**
         * In order to solve the multi dimensional Heston Equation you have to provided
         * some statistical data about the underlying (long-run variance, vol of vol, mean reversion rate). This function allows you to set this data.
         *
         * @param thetas_arg vector of long-run variances, one for each asset.
         * @param kappas_arg vector of mean reversion rates, one for each asset.
         * @param volvols_arg vector of volatility of volatilities, one for each asset.
         * @param rhos correlation matrix. Size needs to be twice the number of assets to handle the two Wiener processes for each asset.
         * @param r market risk-free interest rate scalar.
         */
        virtual void setStochasticData(SGPP::base::DataVector& thetas_arg, SGPP::base::DataVector& kappas_arg, SGPP::base::DataVector& volvols_arg, SGPP::base::DataMatrix& rhos, float_t r);

        void solveImplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveExplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        /**
         * Concrete implementation of the Crank-Nicholson solver. Optional use of implicit euler as a 'kickstart'.
         *
         * @param numTimesteps the number of timesteps to use. The product of this and the timestepsize parameter must be equal to the maturity T.
         * @param timestepsize size of the timestep to take. Assumed constant.
         * @param maxCGIterations the maximum number of iterations to be undertaken by the CG algorithm before aborting.
         * @param epsilonCG residual value below which the CG iterative process stops
         * @param alpha coefficients of the sparse grid basis functions
         * @param NumImEul the number of initial implicit euler steps to take before switching to Crank-Nicholson
         */
        void solveCrankNicolson(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, size_t NumImEul = 0);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType);

        /**
         * Inits the screen object
         */
        virtual void initScreen();

        /**
         * returns the algorithmic dimensions (the dimensions in which the Up Down
         * operations (need for space discretization) should be applied)
         *
         * @return the algorithmic dimensions
         */
        virtual std::vector<size_t> getAlgorithmicDimensions();

        /**
         * sets the algorithmic dimensions (the dimensions in which the Up Down
         * operations (need for space discretization) should be applied)
         *
         * @param newAlgoDims std::vector containing the algorithmic dimensions
         */
        virtual void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

        /**
         *  enables coarsening of grid during solving the Heston
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
        virtual void setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, float_t coarsenThreshold, float_t refineThreshold);

        /**
         * Evaluates the current option value
         * at a point given in Cartesian coordinates
         *
         * @param eval_point the point at with the option price should be determined
         * @param alpha the grid's coefficients
         *
         * @return the option price at the given point
         */
        virtual float_t evalOption(std::vector<float_t>& eval_point, SGPP::base::DataVector& alpha);

        /**
         * This method transforms a point given
         * in Cartesian coordinates into the coordinates used by the
         * current instance of HestonSolver
         *
         * @param point point given in Cartision coordinates that should be transformed
         */
        virtual void transformPoint(SGPP::base::DataVector& point);

        /**
         * Resets the current solving time.
         *
         * use this in order to get the discounting right when using one
         * instance of multiple option pricings
         */
        virtual void resetSolveTime();

        /**
         * gets the number of gridpoints at money
         *
         * Only on Cartesian grids!
         *
         * @param payoffType the payoff type
         * @param strike the option's strike
         * @param eps epsilon to determine the gridpoints, use if at money is not exactly on grid
         * @return number of gridpoints at money
         */
        virtual size_t getGridPointsAtMoney(std::string payoffType, float_t strike, float_t eps = 0.0);

        /**
         * gets the number needed iterations to solve the Heston Equation
         *
         * @return number of iterations needed to solve the Heston Equation, if called before solving 0 is returned
         */
        virtual size_t getNeededIterationsToSolve();

        /**
         * gets needed time in seconds to solve the Heston Equation
         *
         * @return needed time in seconds to solve the Heston Equation, if called before solving 0 is returned
         */
        virtual float_t getNeededTimeToSolve();

        /**
         * gets the number of points in start grid
         *
         * @returns the number of points in start grid, if called before constructing grid, 0 is returned
         */
        virtual size_t getStartInnerGridSize();

        /**
         * gets the number of points in final grid
         *
         * @returns the number of points in final grid, if called before solving, 0 is returned
         */
        virtual size_t getFinalInnerGridSize();

        /**
         * gets the number of average gridpoints
         *
         * @returns the number of average gridpoints, if called before solving, 0 is returned
         */
        virtual size_t getAverageInnerGridSize();

        /**
         * Uses the current grid to evaluate the closed-form Heston surface.
         *
         * @param alpha the vector in which to store the closed-form results
         * @param maturity the option maturity
         */
        void EvaluateHestonExactSurface(SGPP::base::DataVector& alpha, float_t maturity);

        /**
         * Uses the put-call parity to evaluate the exact Heston surface for a put option, based on the closed-form Heston surface for a call option.
         *
         * @param alpha the vector in which to store the results
         * @param maturity the option maturity
         */
        void EvaluateHestonExactSurfacePut(SGPP::base::DataVector& alpha, float_t maturity);

        /**
         * Evaluates the Heston closed-form curve for a constant variance and varying stock price based on the provided grid.
         *
         * @param alpha the vector in which to store the results
         * @param grid1d the grid that specifies the stock price evaluation points
         * @param boundingBox1d bounding box for the grid
         * @param maturity the option maturity
         * @param v the constant variance value for evaluation
         */
        void EvaluateHestonExact1d(SGPP::base::DataVector& alpha, SGPP::base::Grid* grid1d, SGPP::base::BoundingBox* boundingBox1d, float_t maturity, float_t v);

        /**
         * Uses a Black-Scholes solver to evaluate the curve for a constant volatility.
         *
         * @param alpha the vector in which to store the results
         * @param grid1d the grid that specifies the stock price evaluation points
         * @param boundingBox1d bounding box for the grid
         * @param maturity the option maturity
         * @param sigma the constant volatility to use for evaluation
         */
        void EvaluateBsExact1d(SGPP::base::DataVector& alpha, SGPP::base::Grid* grid1d, SGPP::base::BoundingBox* boundingBox1d, float_t maturity, float_t sigma);

        /**
         * Evaluates the closed-form Heston price for a vanilla call.
         *
         * @param S stock price
         * @param v variance
         * @param xi volatility of the volatility
         * @param theta long-run variance
         * @param kappa mean reversion rate
         * @param rho correlation between the stock-price process and the variance process
         * @param r market risk-free interest rate
         * @param T option maturity
         * @param K option strike price
         *
         * @return option price
         */
        float_t EvaluateHestonPriceExact(float_t S, float_t v, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K);

        /*
         * Calls the larger overload (the one with more parameters) of EvaluateHestonPriceExact with the stochastic data already saved in this instance.
         *
         * @param S stock price
         * @param v variance
         */
        float_t EvaluateHestonPriceExact(float_t S, float_t v, float_t maturity);

        /**
         * Evaluates the closed-form Heston price for a vanilla put.
         *
         * @param S stock price
         * @param v variance
         * @param xi volatility of the volatility
         * @param theta long-run variance
         * @param kappa mean reversion rate
         * @param rho correlation between the stock-price process and the variance process
         * @param r market risk-free interest rate
         * @param T option maturity
         * @param K option strike price
         *
         * @return option price
         */
        float_t EvaluateHestonPriceExactPut(float_t S, float_t v, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K);

        /*
         * Calls the larger overload (the one with more parameters) of EvaluateHestonPriceExactPut with the stochastic data already saved in this instance.
         *
         * @param S stock price
         * @param v variance
         */
        float_t EvaluateHestonPriceExactPut(float_t S, float_t v, float_t maturity);

        /**
         * Evaluates the difference surface between the closed-form Heston surface and the closed-form Black-Scholes surface on the current grid.
         *
         * @param alpha vector in which to store the difference results
         * @param maturity the option maturity
         */
        void CompareHestonBsExact(SGPP::base::DataVector& alpha, float_t maturity);

        /**
         * Prints the closed-form BS, Heston and difference curves for a constant variance (and thus constant volatility) to files.
         *
         * @param maturity option maturity
         * @param v variance value (from which volatility is also derived)
         */
        void CompareHestonBs1d(float_t maturity, float_t v);

        /**
         * Perform a single step of the Gauss-Lobatto integration for the Heston closed-form integral functions.
         * Adapted from the implementation found at http://www.bnikolic.co.uk/nqm/1dinteg/gausslobatto.html
         *
         * @param a lower integration limit
         * @param b upper integration limit
         * @param fa value of the function at the lower limit (used to save an evaluation when refinement is used)
         * @param fb value of the function at the upper limit (used to save an evaluation when refinement is used)
         * @param neval number of evaluations made so far
         * @param maxeval maximum number of evaluations which should not be exceeded
         * @param acc required accuracy expressed in units of std::numeric_limits<float_t>::epsilon(). This allows less-than comparison by using addition and equality.
         * @param xi Heston parameter
         * @param theta Heston parameter
         * @param kappa Heston parameter
         * @param rho Heston parameter
         * @param r discount rate
         * @param T Heston parameter
         * @param K Heston parameter
         * @param S Heston parameter
         * @param v Heston parameter
         * @param type Heston parameter
         *
         * @return integral step value
         */
        float_t GaussLobattoIntStep(float_t a, float_t b, float_t fa, float_t fb, size_t& neval, size_t maxeval, float_t acc, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K, float_t S, float_t v, int type);

        /**
         * Perform Gauss-Lobatto integration for the Heston closed-form integral functions.
         * Adapted from the implementation found at http://www.bnikolic.co.uk/nqm/1dinteg/gausslobatto.html
         *
         * @param a the lower integration limit
         * @param b the upper integration limit
         * @param abstol absolute tolerance - integration stops when the error estimate is smaller than this
         * @param maxeval maximum of evaluations to make. If this number of evaluations is made without reaching the required accuracy, an exception of type std::runtime_error is thrown.
         * @param xi Heston parameter
         * @param theta Heston parameter
         * @param kappa Heston parameter
         * @param rho Heston parameter
         * @param r discount rate
         * @param T Heston parameter
         * @param K Heston parameter
         * @param S Heston parameter
         * @param v Heston parameter
         * @param type Heston parameter
         *
         * @return integral value
         */
        float_t GaussLobattoInt(float_t a, float_t b, float_t abstol, size_t maxeval, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K, float_t S, float_t v, int type);


        /**
         * Compares a numerical Heston solution to the exact solution and writes the difference to a file.
         *
         * @param solution numerical Heston solution
         * @param exact (closed-form) Heston solution
         * @param filename name of the file in which to write the difference
         * @param PointsPerDimension granularity of the difference measurements
         */
        void CompareHestonSolutionToExact(SGPP::base::DataVector* solution, SGPP::base::DataVector* exact, std::string filename, size_t PointsPerDimension);

        /**
         * Uses an evaluation operator to get the interpolated option value for a particular s and v based on the given numerical solution.
         *
         * @param s stock-price for the evaluation point
         * @param v variance for the evaluation point
         * @param alphaVec the numerical solution on which to do the evaluation
         *
         * @return the evaluated option price
         */
        float_t EvalSinglePoint1Asset(float_t s, float_t v, SGPP::base::DataVector& alphaVec);

        /**
         * Gets the closed-form Black-Scholes surface on the current grid. Uses the equivalent volatility of the variance values. Vanilla option only.
         *
         * @param alphaBS vector in which to store the exact BS solution values
         * @param maturity the option maturity
         */
        void GetBsExactSolution(SGPP::base::DataVector& alphaBS, float_t maturity);

        /**
         * Computes the difference between the a provided numerical Heston solution and the equivalent Black-Scholes solution on the current grid.
         *
         * @param alphaHestonNumeric numeric Heston solution
         * @param alphaBS vector in which to store the equivalent BS solution
         * @param error vector in which to store the difference values
         * @param maturity the option maturity
         */
        void CompareHestonNumericToBsExact(SGPP::base::DataVector& alphaHestonNumeric, SGPP::base::DataVector& alphaBS, SGPP::base::DataVector& error, float_t maturity);
    };

  }
}

#endif /* HESTONSOLVER_HPP */