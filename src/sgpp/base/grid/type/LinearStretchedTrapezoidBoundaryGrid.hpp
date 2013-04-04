/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP
#define LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg {
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

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP */
