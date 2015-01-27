// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARGENERALIZEDTRUNCATEDBOUNDARYGRID_HPP_
#define LINEARGENERALIZEDTRUNCATEDBOUNDARYGRID_HPP_
#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions with boundaries, pentagon cut
     * Generalization of the LinearBoundary and LinearTruncatedBoundary Grids
     * The sparse grid does contain all fullgrids with |l|<a given level, and l_min>l_user
     * For l_user=0 we obtain the LinearBoundaryGrid and for l_user=1 we obtain the linear truncated boundary grid
     */
    class LinearGeneralizedTruncatedBoundaryGrid : public Grid {
      protected:
        LinearGeneralizedTruncatedBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Truncated Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        LinearGeneralizedTruncatedBoundaryGrid(size_t dim);

        /**
         * Constructor Linear Truncated Boundary Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearGeneralizedTruncatedBoundaryGrid(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~LinearGeneralizedTruncatedBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARGENERALIZEDTRUNCATEDBOUNDARYGRID_HPP_ */
