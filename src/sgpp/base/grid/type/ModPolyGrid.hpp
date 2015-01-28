/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef MODPOLYGRID_HPP
#define MODPOLYGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg {
  namespace base {

    /**
     * grid with modified polynomial base functions
     */
    class ModPolyGrid : public Grid {
      protected:
        ModPolyGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified polynomial base functions
         *
         * @param dim the dimension of the grid
         * @param degree the max. polynom's degree
         */
        ModPolyGrid(size_t dim, size_t degree);

        /**
         * Destructor
         */
        virtual ~ModPolyGrid();

        virtual const char* getType();
        virtual void serialize(std::ostream& ostr);

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
        virtual size_t getDegree() const;

      protected:
        /// max. polynom's degree
        size_t degree;
    };

  }
}

#endif /* MODPOLYGRID_HPP */
