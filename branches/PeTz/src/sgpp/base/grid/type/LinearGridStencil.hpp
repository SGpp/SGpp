/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Gerrit Buse (buse@in.tum.de)

#ifndef LINEARGRIDSTENCIL_HPP
#define LINEARGRIDSTENCIL_HPP

#include "base/grid/type/GridStencil.hpp"
#include "base/grid/common/BoundingBox.hpp"

#include <iostream>

namespace sg {
  namespace base {

    /**
     * grid with linear base functions
     */
    class LinearGridStencil : public GridStencil {
      protected:
        LinearGridStencil(std::istream& istr);

      public:
        /**
         * Constructor Linear Grid without boundaries
         *
         * @param dim the dimension of the grid
         */
        LinearGridStencil(size_t dim);

        /**
         * Constructor Linear Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearGridStencil(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~LinearGridStencil();

        virtual const char* getType();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);


    };

  }
}

#endif /* LINEARGRIDSTENCIL_HPP */
