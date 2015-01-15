/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef MODLINEARGRID_HPP
#define MODLINEARGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg {
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

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

    };

  }
}

#endif /* MODLINEARGRID_HPP */
