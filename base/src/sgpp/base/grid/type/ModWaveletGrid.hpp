// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODWAVELETGRID_HPP
#define MODWAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with modified wavelet base functions
     */
    class ModWaveletGrid : public Grid {
      protected:
        ModWaveletGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified wavelet base functions
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
