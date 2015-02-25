// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PERIODICGRID_HPP
#define PERIODICGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with modified linear base functions
     */
    class PeriodicGrid : public Grid {
      protected:
        PeriodicGrid(std::istream& istr);

      public:
        /**
         * Constructor modified linear grid
         *
         * @param dim the dimension of the grid
         */
        PeriodicGrid(size_t dim);

        /**
         * Destructor
         */
        virtual ~PeriodicGrid();

        virtual const char* getType();

        virtual GridGenerator* createGridGenerator();

        virtual const SBasis& getBasis();

        static Grid* unserialize(std::istream& istr);

    };

  }
}

#endif /* PERIODICGRID_HPP */