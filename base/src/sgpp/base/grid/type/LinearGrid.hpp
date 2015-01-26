// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LINEARGRID_HPP
#define LINEARGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions
     */
    class LinearGrid : public Grid {
      protected:
        LinearGrid(std::istream& istr);

      public:
        /**
         * Constructor Linear Grid without boundaries
         *
         * @param dim the dimension of the grid
         */
        LinearGrid(size_t dim);

        /**
         * Constructor Linear Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearGrid(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~LinearGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARGRID_HPP */