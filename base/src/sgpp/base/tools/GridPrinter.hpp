/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#ifndef GRIDPRINTER_HPP
#define GRIDPRINTER_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements a utility that allows you to print a grid
     * to file. These files can be plotted with gnuplot.
     */
    class GridPrinter {
      protected:
        /// Pointer to the grid Object, which should be printed
        Grid* myGrid;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid Reference to a Spare Grid, that should be printed
         */
        GridPrinter(Grid& SparseGrid);

        /**
         * Std-Destructor
         */
        virtual ~GridPrinter();

        /**
         * Print the grid points in level/index format to a file; front end
         *
         * @param tFilename absolute path to the file the grid is written into
         */
        virtual void printLevelIndexGrid(std::string tFilename);

        /**
         * Print the grid with its function to a file; front end
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absolute path to the file the grid is written into
         * @param PointsPerDimension specifies how many functions evaluations in every dimension should be calculated
         */
        virtual void printGrid(DataVector& alpha, std::string tFilename,
                               size_t PointsPerDimension);

        /**
         * Print the grid with its function to a file; front end
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absolute path to the file the grid is written into
         * @param GridArea The area in which the function should be plotted
         * @param PointsPerDimension specifies how many functions evaluations in every dimension should be calculated
         */
        virtual void printGridDomain(DataVector& alpha, std::string tFilename,
                                     BoundingBox& GridArea, size_t PointsPerDimension);

        /**
         * Prints the Grid Points of the Sparse Grid either with their node basis value
         * or their hierarchical surplus
         *
         * @param alpha the coefficients of the grid's ansatzfunctions
         * @param tFilename absoulte path to the file the grid is written into
         * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
         */
        virtual void printSparseGrid(DataVector& alpha, std::string tFilename,
                                     bool bSurplus);

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
        virtual void printSparseGridExpTransform(DataVector& alpha,
            std::string tFilename, bool bSurplus);
    };

  }
}

#endif /* GRIDPRINTER */
