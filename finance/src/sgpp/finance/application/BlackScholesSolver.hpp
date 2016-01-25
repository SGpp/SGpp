// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESSOLVER_HPP
#define BLACKSCHOLESSOLVER_HPP


#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/solver/ODESolver.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>

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
#include "../../../../../base/src/sgpp/base/grid/type/LinearBoundaryGrid.hpp"


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
    class BlackScholesSolver : public SGPP::pde::ParabolicPDESolver {
      protected:
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
        /// identifies if the Black Scholes Equation should be solved by using a principal axis transformation
        bool usePAT;
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
        /// Type of the Option to solve
        std::string tBoundaryType;
        /// Eigenvectors of the co-variance matrix
        SGPP::base::DataMatrix* eigvec_covar;
        /// Eigenvalues of the co-variance matrix
        SGPP::base::DataVector* eigval_covar;
        /// mu hat, tanslation coefficient needed if PAT is used
        SGPP::base::DataVector* mu_hat;
        /// stores the current time until which the option has been solved
        float_t current_time;
        /// stores the strike of the current option
        float_t dStrike;
        /// stores the option type of the current option
        std::string payoffType;

        /**
         * returns the option value (payoff value) for an European call option
         *
         * @param assetValue the current asset's value
         * @param strike the strike price of the option
         *
         * @return the call premium
         */
        virtual float_t get1DEuroCallPayoffValue(float_t assetValue, float_t strike);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option.
         * The grid is initialized based on Cartesian coordinates!
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initCartesianGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         * The grid is initialized based on log-transformed coordinates!
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initLogTransformedGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         * The grid is initialized based on log-transformed and a principal axis transformation coordinates!
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initPATTransformedGridWithPayoff(SGPP::base::DataVector& alpha, float_t strike, std::string payoffType);

        /**
         * This function calculates for every grid point the value
         * of a normal distribution given by norm_mu and norm_sigma.
         * The result is stored dehierarchized in alpha.
         *
         * This method is overwritten in order to support grids with logarithmic coordinates.
         *
         * @param alpha contains dehierarchized sparse grid coefficients containing the values of the multi dimensional normal distribution after call
         * @param norm_mu the expected values of the normal distribution for every grid dimension
         * @param norm_sigma the standard deviation of the normal distribution for every grid dimension
         */
        virtual void getGridNormalDistribution(SGPP::base::DataVector& alpha, std::vector<float_t>& norm_mu, std::vector<float_t>& norm_sigma);

      public:
        /**
         * Std-Constructor of the solver
         *
         * @param useLogTransform speciefies if a log transformed formulation should be used for solving BlackScholes Equation
         * @param usePAT speciefies if a principal axis transformation (also enabling a log-transformation) should be used for solving BlackScholes Equation
         */
        BlackScholesSolver(bool useLogTransform = false, bool usePAT = false);

        /**
         * Std-Destructor of the solver
         */
        virtual ~BlackScholesSolver();

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
         * In order to solve the multi dimensional Black Scholes Equation you have to provided
         * some statistical data about the underlying (assets' weight, standard deviation
         * and the correlation between them). This function allows you to set this data.
         *
         * @param mus a SGPP::base::DataVector that contains the underlyings' weight
         * @param sigmas a SGPP::base::DataVector that contains the underlyings' standard deviations
         * @param rhos a SGPP::base::DataMatrix that contains the correlations between the underlyings
         * @param r the riskfree rate used in the market model
         */
        virtual void setStochasticData(SGPP::base::DataVector& mus, SGPP::base::DataVector& sigmas, SGPP::base::DataMatrix& rhos, float_t r);

        void solveImplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveExplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

        void solveCrankNicolson(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, size_t NumImEul = 0);

        void solveX(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false, void* myODESolverV = NULL, std::string Solver = "ImEul");

        void solveSC(std::string Solver, size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        void solveAdamsBashforth(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        void solveSCAC(size_t numTimesteps, float_t timestepsize, float_t epsilon, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        void solveSCH(size_t numTimesteps, float_t timestepsize, float_t epsilon, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        void solveSCBDF(size_t numTimesteps, float_t timestepsize, float_t epsilon, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        void solveSCEJ(size_t numTimesteps, float_t timestepsize, float_t epsilon, float_t myAlpha, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        //void solveSCMC(size_t numTimesteps, float_t timestepsize, float_t epsilon, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

        /**
         * evaluates the analytic solution of the 1D Black Scholes equation
         * @param stock  stock value (point of evaluation)
         * @param isCall set this to true to calculate call, false calculates put
         * @param t time to maturity
         * @param vola volatility
         * @param r risk-free rate
         * @param strike the strike price of the option
         *
         * @returns the option price for the given stock value
         */
        float_t getAnalyticSolution1D(float_t stock, bool isCall, float_t t, float_t vola, float_t r, float_t strike);

        /**
         * Solves the closed form of the Black Scholes equation, the Black Scholes
         * formular. It evaluates the Black Scholes formular in a Stock Price Range
         * from 0.0 to maxStock, by increasing the stock price in every step by
         * a given (small) values, so the analytical solution of the PDE can
         * be determined and compared.
         *
         * @param premiums the result vector, here the combinations of stock price and premium are stored
         * @param minStock the minimum stock regarded in these calculations
         * @param maxStock the maximum stock regarded in these calculations
         * @param StockInc the increase of the stockprice in one step
         * @param strike the strike price of the Option
         * @param t time to maturity
         * @param isCall set this to true to calculate call, false calculates put
         */
        void solve1DAnalytic(std::vector< std::pair<float_t, float_t> >& premiums, float_t minStock, float_t maxStock, float_t StockInc, float_t strike, float_t t, bool isCall);

        /**
         * Writes the premiums into a file that can be easily plot with gnuplot
         *
         * @param premiums the result vector, here the combinations of stock price and premium are stored
         * @param tfilename absolute path to file into which the grid's evaluation is written
         */
        void print1DAnalytic(std::vector< std::pair<float_t, float_t> >& premiums, std::string tfilename);

        /**
         *  computes the relative error between the solution and the exact analytic solution for the 1-dimensional Black-Schoesl equation
         *
         *  @param alpha_analytic data vector with the analytic solution
         *  @param strike strike price of the option
         *  @param t maturity time
         *  @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         *  @param hierarchized flag whether values should be hierarchized (true=hierarchized, false=dehierarchized)
         */
        virtual void getAnalyticAlpha1D(SGPP::base::DataVector& alpha_analytic, float_t strike, float_t t, std::string payoffType, bool hierarchized);

        /**
         * Evaluates the analytic solution of the 1d Black Scholes equation
         * on different points specified in EvaluationPoints and stores the result into FunctionValues.
         *
         * @param AnalyticOptionPrices SGPP::base::DataVector into the which the result of function's evaluation is stored
         * @param EvaluationPoints SGPP::base::DataMatrix that contains the points at which the sparse grid's function is evaluated
         * @param strike the strike price of the Option
         * @param vola volatility
         * @param r risk-free rate
         * @param t time to maturity
         * @param isCall set this to true to calculate call, false calculates put
         */
        void evaluate1DAnalyticCuboid(SGPP::base::DataVector& AnalyticOptionPrices, SGPP::base::DataMatrix& EvaluationPoints, float_t strike, float_t vola, float_t r, float_t t, bool isCall);

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
        virtual void setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, float_t coarsenThreshold, float_t refineThreshold);

        /**
         * prints the 2D interpolation error at money into a file. This file is plotable via gnuplot. A bounding
         * box [0,x] X [0,y] is assumed.
         *
         * Only on Cartesian grids!
         *
         * @param alpha the sparse grid's coefficients
         * @param tFilename the name of file contain the interpolation error
         * @param numTestpoints Number of equal distribute testpoints at money
         * @param strike the option's strike
         */
        virtual void printPayoffInterpolationError2D(SGPP::base::DataVector& alpha, std::string tFilename, size_t numTestpoints, float_t strike);

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
         * current instance of BlackScholesSolver
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
         * Prints the SGPP::base::Grid Points of the Sparse SGPP::base::Grid either with their node basis value
         * or their hierarchical surplus
         *
         * This function is available for all dimensions
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tfilename absoulte path to the file the grid is written into
         * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
         */
        void printSparseGridPAT(SGPP::base::DataVector& alpha, std::string tfilename, bool bSurplus) const;

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
         * gets the number needed iterations to solve Black Scholes Equation
         *
         * @return number of iterations needed to solve Black Scholes Equation, if called before solving 0 is returned
         */
        virtual size_t getNeededIterationsToSolve();

        /**
         * gets needed time in seconds to solve Black Scholes Equation
         *
         * @return needed time in seconds to solve Black Scholes Equation, if called before solving 0 is returned
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
         * Routine to export the RHS of the inner system which has to be
         * solved in order to solve the Poisson equation
         *
         * @param alpha the start solution
         * @param tFilename file into which the rhs is written
         * @param timestepsize the size of the timesteps
         */
        void storeInnerRHS(SGPP::base::DataVector& alpha, std::string tFilename, float_t timestepsize);

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
        void storeInnerSolution(SGPP::base::DataVector& alpha, size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, std::string tFilename);
    };

  }
}

#endif /* BLACKSCHOLESSOLVER_HPP */