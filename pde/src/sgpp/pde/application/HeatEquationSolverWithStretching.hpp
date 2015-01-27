// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef HEATEQUATIONSOLVERWITHSTRETCHING_HPP
#define HEATEQUATIONSOLVERWITHSTRETCHING_HPP


#include <sgpp/pde/application/ParabolicPDESolver.hpp>

//#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
//#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/base/tools/GridPrinterForStretching.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include <sgpp/globaldef.hpp>
#include "../../../../../base/src/sgpp/base/grid/type/LinearStretchedTruncatedBoundaryGrid.hpp"


namespace SGPP {
  namespace pde {

    /**
     * This class provides a simple-to-use solver of the multi dimensional
     * Heat Equation that uses Sparse Grids.
     *
     * The class's aim is, to hide all complex details of solving the
     * Heat Equation with Sparse Grids!
     *
     * @version $HEAD$
     */
    class HeatEquationSolverWithStretching : public ParabolicPDESolver {
      private:
        /// the heat coefficient
        double a;
        /// screen object used in this solver
        SGPP::base::ScreenOutput* myScreen;
        ////BoundingBox replacement
        SGPP::base::Stretching* myStretching;

      public:
        /**
         * Std-Constructor of the solver
         */
        HeatEquationSolverWithStretching();

        /**
         * Std-Destructor of the solver
         */
        virtual ~HeatEquationSolverWithStretching();

        void constructGrid(SGPP::base::Stretching& myStretching, int level);

        void constructGrid(SGPP::base::BoundingBox& myStretching, int level);

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
         * Inits the grid in the middle of the whole domain with one single heat
         *
         * alpha reference to the coefficients vector
         * heat the value of the heat in the middle of the domain
         */
        //  void initGridWithSingleHeat(SGPP::base::DataVector& alpha, double heat);

        /**
         * Inits the grid in the middle the domain with an smooth heat distribution that the
         * normal distribution formula
         *
         * @param alpha reference to the coefficients vector
         * @param mu the exspected value of the normal distribution
         * @param sigma the sigma of the normal distribution
         * @param factor a factor that is used to stretch the function values
         */
        void initGridWithSmoothHeat(SGPP::base::DataVector& alpha, double mu, double sigma, double factor);

        /**
         * Inits the grid with a constant heat
         *
         * alpha reference to the coefficients vector
         * constHeat the temperature of the constant heat
         */
        //  void initGridWithConstantHeat(SGPP::base::DataVector& alpha, double constHeat);

        /**
         * Inits the screen object
         */
        void initScreen();

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
        void printGrid(SGPP::base::DataVector& alpha, double PointesPerDimension, std::string tfilename) const;

        /**
         This function is a placeholder, is not used.
         */
        void printGridDomain(SGPP::base::DataVector& alpha, double PointesPerDimension, SGPP::base::BoundingBox& GridArea, std::string tfilename) const;

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
        void printGridDomainStretching(SGPP::base::DataVector& alpha, double PointesPerDimension, SGPP::base::Stretching& GridArea, std::string tfilename) const;

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
        void printSparseGrid(SGPP::base::DataVector& alpha, std::string tfilename, bool bSurplus) const;

        /**
         * Prints the SGPP::base::Grid Points of the Sparse SGPP::base::Grid either with their node basis value
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
        void printSparseGridExpTransform(SGPP::base::DataVector& alpha, std::string tfilename, bool bSurplus) const;

    };

  }
}

#endif /* HEATEQUATIONSOLVERWITHSTRETCHING_HPP */