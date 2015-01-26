/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef BLACKSCHOLESSOLVERWITHSTRETCHING_HPP
#define BLACKSCHOLESSOLVERWITHSTRETCHING_HPP

#include "finance/application/BlackScholesSolver.hpp"
#include "base/grid/type/LinearStretchedTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/LinearStretchedGrid.hpp"
#include "base/grid/common/Stretching.hpp"
#include "base/tools/GridPrinterForStretching.hpp"

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
    class BlackScholesSolverWithStretching : public BlackScholesSolver {
      private:
        /// Stores Pointer to the sg::base::Grid's sg::base::Stretching
        sg::base::Stretching* myStretching;

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option.
         * The grid is initialized based on Cartesian coordinates!
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initCartesianGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         * The grid is initialized based on log-transformed coordinates!
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initLogTransformedGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType);

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
        virtual void getGridNormalDistribution(sg::base::DataVector& alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma);

      public:
        /**
         * Std-Constructor of the solver
         *
         * @param useLogTransform specifies if a log transformed formulation should be used for solving BlackScholes Equation
         * @param OptionType possible values "all" and "European", if "European" is choose a solver with fix Dirichlet boundaries is selected
         */
        BlackScholesSolverWithStretching(bool useLogTransform = false, std::string OptionType = "all");

        /**
         * Std-Destructor of the solver
         */
        ~BlackScholesSolverWithStretching();

        void constructGridStretching(sg::base::Stretching& myStretching, int level);

        void constructGrid(sg::base::BoundingBox& myBoundingBox, size_t level);

        /**
         * This function tries to refine the grid such that
         * most of the grid points are used for interpolation of the singularity. So this grid
         * is able to approximate the start solution better.
         *
         * After refining the grid the payoff function is applied to the grid.
         *
         * Only on Cartesian grids!
         *
         * @param alpha reference to a sg::base::DataVector object that contains the gird ansatzfunction's coefficients
         * @param strike containing the option's strike
         * @param payoffType the type of payoff Function used ONLY supported: avgM
         * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
         */
        virtual void refineInitialGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance);

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
         * @param alpha reference to a sg::base::DataVector object that contains the gird ansatzfunction's coefficients
         * @param strike containing the option's strike
         * @param payoffType the type of payoff Function used ONLY supported: avgM
         * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
         * @param maxLevel maximum level of refinement
         */
        virtual void refineInitialGridWithPayoffToMaxLevel(sg::base::DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance, sg::base::GridIndex::level_type maxLevel);

        /**
         * Inits the alpha vector with a payoff function of an European call option or put option
         *
         * @param alpha the coefficient vector of the grid's ansatzfunctions
         * @param strike the option's strike
         * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         */
        virtual void initGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType);

        /**
         *  computes the relative error between the solution and the exact analytic solution for the 1-dimensional Black-Schoesl equation
         *
         *  @param alpha_analytic data vector with the analytic solution
         *  @param strike strike price of the option
         *  @param t maturity time
         *  @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
         *  @param hierarchized flag whether values should be hierarchized (true=hierarchized, false=dehierarchized)
         */
        virtual void getAnalyticAlpha1D(sg::base::DataVector& alpha_analytic, double strike, double t, std::string payoffType, bool hierarchized);

        /**
         * Inits the screen object
         */
        virtual void initScreen();


        /**
         * prints the 2D interpolation error at money into a file. This file is plotable via gnuplot. A bounding
         * box [0,x] X [0,y] is assumed.
         *
         * Only on Cartesian grids!
         *
         * @param alpha the sparse grid's coefficients
         * @param tFilename the name of file contain the interpolation error
         * @param numTestpoints Number of equal distribute testpoints at the money
         * @param strike the option's strike
         */
        virtual void printPayoffInterpolationError2D(sg::base::DataVector& alpha, std::string tFilename, size_t numTestpoints, double strike);

        /**
         * gets the number of gridpoints at money
         *
         * Only on Cartesian grids!
         *
         * @param payoffType the payoff type
         * @param strike the option's strike
         * @param eps epsilon to determine the gridpoints, use if at money is not exactly on grid
         */
        virtual size_t getGridPointsAtMoney(std::string payoffType, double strike, double eps = 0.0);


        /**
         * This is some kind of debug functionality. It writes a file,
         * that can be used with gnuplot the print the grid.
         *
         * Is only implemented for 1D and 2D grids!
         *
         * @param alpha the coefficients of the Sparse Gird's basis functions
         * @param PointesPerDimension the distance between evaluation points
         * @param tfilename absolute path to file into which the grid's evaluation is written
         */
        virtual void printGrid(sg::base::DataVector& alpha, size_t PointesPerDimension, std::string tfilename) const;

        /**
         * This is not used, throws exception to inform about the function printGridDomainStretching
         */
        virtual void printGridDomain(sg::base::DataVector& alpha, size_t PointesPerDimension, sg::base::BoundingBox& GridArea, std::string tfilename) const;

        /**
         * This is some kind of debug functionality. It writes a file,
         * that can be used with gnuplot the print the grid.
         *
         * Is only implemented for 2D grids!
         *
         * @param alpha the coefficients of the Sparse Gird's basis functions
         * @param PointesPerDimension the distance between evaluation points
         * @param GridArea the area in which the function should be plotted
         * @param tfilename absolute path to file into which the grid's evaluation is written
         */
        virtual void printGridDomainStretching(sg::base::DataVector& alpha, size_t PointesPerDimension, sg::base::Stretching& GridArea, std::string tfilename) const;

        /**
         * Prints the sg::base::Grid Points of the Sparse sg::base::Grid either with their node basis value
         * or their hierarchical surplus
         *
         * This function is available for all dimensions
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tfilename absoulte path to the file the grid is written into
         * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
         */
        virtual void printSparseGrid(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const;

        /**
         * Prints the sg::base::Grid Points of the Sparse sg::base::Grid either with their node basis value
         * or their hierarchical surplus
         *
         * This function is available for all dimensions.
         *
         * The coordinates of the grid points are pushed the exp function. So
         * log transformed grids can be plotted in cartesion coordinates.
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tfilename absoulte path to the file the grid is written into
         * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
         */
        virtual void printSparseGridExpTransform(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const;


    };

  }
}
#endif /* BLACKSCHOLESSOLVERWITHSTRETCHING_HPP */
