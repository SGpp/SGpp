// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODLINEARGRID_HPP
#define MODLINEARGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with modified linear base functions
     */
    class ModLinearGrid : public Grid {
      protected:
        ModLinearGrid(std::istream& istr);

      public:
        /**
         * Constructor modified linear grid
         *
         * @param dim the dimension of the grid
         */
        ModLinearGrid(size_t dim);

        /**
         * Destructor
         */
        virtual ~ModLinearGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

    };

  }
}

#endif /* MODLINEARGRID_HPP */