// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARSTRETCHEDTRUNCATEDBOUNDARYGRID_HPP
#define LINEARSTRETCHEDTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions with boundaries, pentagon cut
     */
    class LinearStretchedTruncatedBoundaryGrid : public Grid {
      protected:
        LinearStretchedTruncatedBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Truncated Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        LinearStretchedTruncatedBoundaryGrid(size_t dim);

        /**
         * Constructor Linear Truncated Boundary Grid
         *
         * @param BB the Stretching of the grid
         */
        LinearStretchedTruncatedBoundaryGrid(Stretching& BB);

        /**
         * Destructor
         */
        virtual ~LinearStretchedTruncatedBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARSTRETCHEDTRUNCATEDBOUNDARYGRID_HPP */
