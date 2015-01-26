/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Gerrit Buse (buse@in.tum.de)

#ifndef MODLINEARGRIDSTENCIL_HPP
#define MODLINEARGRIDSTENCIL_HPP

#include <sgpp/base/grid/type/GridStencil.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions
     */
    class ModLinearGridStencil : public GridStencil {
      protected:
        ModLinearGridStencil(std::istream& istr);

      public:
        /**
         * Constructor Linear Grid without boundaries
         *
         * @param dim the dimension of the grid
         */
        ModLinearGridStencil(size_t dim);

        /**
         * Constructor Linear Grid
         *
         * @param BB the BoundingBox of the grid
         */
        ModLinearGridStencil(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~ModLinearGridStencil();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);


    };

  }
}

#endif /* MODLINEARGRIDSTENCIL_HPP */
