/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef MODWAVELETGRID_HPP
#define MODWAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

namespace sg {
  namespace base {

    /**
     * grid with modified polynomial base functions
     */
    class ModWaveletGrid : public Grid {
      protected:
        ModWaveletGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified polynomial base functions
         *
         * @param dim the dimension of the grid
         */
        ModWaveletGrid(size_t dim);

        /**
         * Destructor
         */
        virtual ~ModWaveletGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

    };

  }
}

#endif /* MODWAVELETGRID_HPP */
