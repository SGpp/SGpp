// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARBOUNDARYGRID_HPP
#define LINEARBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions with boundaries
     */
    class LinearL0BoundaryGrid : public Grid {
      protected:
        LinearL0BoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor for the Linear Boundary Grid
         *
         * @param dim the dimension of the grid
         */
        LinearL0BoundaryGrid(size_t dim);

        /**
         * Destructor
         */
        virtual ~LinearL0BoundaryGrid();

        virtual SGPP::base::GridType getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARBOUNDARYGRID_HPP */
