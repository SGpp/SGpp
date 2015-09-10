// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARCLENSHAWCURTISGRID_HPP
#define LINEARCLENSHAWCURTISGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with Clenshaw-Curtis linear base functions with boundaries, pentagon cut
     */
    class LinearClenshawCurtisGrid : public Grid {
      protected:
        LinearClenshawCurtisGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Truncated Boundary Clenshaw-Curtis Grid
         *
         * @param dim the dimension of the grid
         */
        LinearClenshawCurtisGrid(size_t dim);

        /**
         * Constructor Linear Truncated Boundary Clenshaw-Curtis Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearClenshawCurtisGrid(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~LinearClenshawCurtisGrid();

        virtual SGPP::base::GridType getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARCLENSHAWCURTISGRID_HPP */
