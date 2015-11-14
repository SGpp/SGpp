// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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

        virtual SGPP::base::GridType getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);


    };

  }
}

#endif /* MODLINEARGRIDSTENCIL_HPP */
