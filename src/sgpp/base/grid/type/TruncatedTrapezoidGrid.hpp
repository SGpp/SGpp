/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author aliz eva

#ifndef TRUNCATEDTRAPEZOIDGRID_HPP_
#define TRUNCATEDTRAPEZOIDGRID_HPP_
#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg {
  namespace base {

    /**
     * grid with linear base functions with boundaries, pentagon cut
     * Generalization of the LinearBoundary and LinearTrapezoidBoundary Grids
     * The sparse grid does contain all fullgrids with |l|<a given level, and l_min>l_user
     * For l_user=0 we obtain the LinearBoundaryGrid and for l_user=1 we obtain the linear trapezoid boundary grid
     */
    class TruncatedTrapezoidGrid : public Grid {
      protected:
        TruncatedTrapezoidGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Trapezoid Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        TruncatedTrapezoidGrid(size_t dim);

        /**
         * Constructor Linear Trapezoid Boundary Grid
         *
         * @param BB the BoundingBox of the grid
         */
        TruncatedTrapezoidGrid(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~TruncatedTrapezoidGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* TRUNCATEDTRAPEZOIDGRID_HPP_ */
