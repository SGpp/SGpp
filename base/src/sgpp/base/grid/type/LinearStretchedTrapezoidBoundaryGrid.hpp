// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP
#define LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions with boundaries, pentagon cut
     */
    class LinearStretchedTrapezoidBoundaryGrid : public Grid {
      protected:
        LinearStretchedTrapezoidBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Trapezoid Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        LinearStretchedTrapezoidBoundaryGrid(size_t dim);

        /**
         * Constructor Linear Trapezoid Boundary Grid
         *
         * @param BB the Stretching of the grid
         */
        LinearStretchedTrapezoidBoundaryGrid(Stretching& BB);

        /**
         * Destructor
         */
        virtual ~LinearStretchedTrapezoidBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP */