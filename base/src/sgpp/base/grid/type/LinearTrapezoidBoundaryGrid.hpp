// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LINEARTRAPEZOIDBOUNDARYGRID_HPP
#define LINEARTRAPEZOIDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions with boundaries, pentagon cut
     */
    class LinearTrapezoidBoundaryGrid : public Grid {
      protected:
        LinearTrapezoidBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Trapezoid Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        LinearTrapezoidBoundaryGrid(size_t dim);

        /**
         * Constructor Linear Trapezoid Boundary Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearTrapezoidBoundaryGrid(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~LinearTrapezoidBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARTRAPEZOIDBOUNDARYGRID_HPP */