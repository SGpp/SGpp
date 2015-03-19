// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARTRUNCATEDBOUNDARYGRID_HPP
#define LINEARTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions with boundaries, pentagon cut
     */
    class LinearTruncatedBoundaryGrid : public Grid {
      protected:
        LinearTruncatedBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Truncated Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        LinearTruncatedBoundaryGrid(size_t dim);

        /**
         * Constructor Linear Truncated Boundary Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearTruncatedBoundaryGrid(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~LinearTruncatedBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARTRUNCATEDBOUNDARYGRID_HPP */
