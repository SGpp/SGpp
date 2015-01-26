/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef PERIODICGRID_HPP
#define PERIODICGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

namespace sg {
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
