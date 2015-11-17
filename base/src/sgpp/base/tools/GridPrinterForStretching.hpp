// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRIDPRINTERFORSTRETCHING_HPP
#define GRIDPRINTERFORSTRETCHING_HPP

#include <sgpp/base/tools/GridPrinter.hpp>

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements a utility that allows you to print a grid
     * to file. These files can be plotted with gnuplot.
     */
    class GridPrinterForStretching: public GridPrinter {
      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid Reference to a Spare Grid, that should be printed
         */
        GridPrinterForStretching(Grid& SparseGrid);

        /**
         * Std-Destructor
         */
        virtual ~GridPrinterForStretching() override;

        /**
         * Print the grid with its function to a file; front end
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absolute path to the file the grid is written into
         * @param PointsPerDimension specifies how many functions evaluations in every dimension should be calculated
         */
        virtual void printGrid(DataVector& alpha, std::string tFilename, size_t PointsPerDimension) override;


        /**
         * This function is not used for stretching grid printing, use printGridDomainStretching instead

         */
        virtual void printGridDomain(DataVector& alpha, std::string tFilename, BoundingBox& GridArea, size_t PointsPerDimension) override;
        /**
         * Print the grid with its function to a file; front end
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absolute path to the file the grid is written into
         * @param GridArea The area in which the function should be plotted
         * @param PointsPerDimension specifies how many functions evaluations in every dimension should be calculated
         */
        void printGridDomainStretching(DataVector& alpha, std::string tFilename, Stretching& GridArea, size_t PointsPerDimension);

        /**
         * Prints the Grid Points of the Sparse Grid either with their node basis value
         * or their hierarchical surplus
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absoulte path to the file the grid is written into
         * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
         */
        virtual void printSparseGrid(DataVector& alpha, std::string tFilename, bool bSurplus) override;

        /**
         * Prints the Grid Points of the Sparse Grid either with their node basis value
         * or their hierarchical surplus.
         *
         * The coordinates of the grid points are pushed the exp function. So
         * log transformed grids can be plotted in cartesion coordinates.
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absoulte path to the file the grid is written into
         * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
         */
        virtual void printSparseGridExpTransform(DataVector& alpha, std::string tFilename, bool bSurplus) override;
    };

  }
}

#endif /* GRIDPRINTERFORSTRETCHING */
