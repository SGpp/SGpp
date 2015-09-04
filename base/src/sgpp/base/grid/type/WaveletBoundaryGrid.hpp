// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELETTRUNCATEDBOUNDARYGRID_HPP
#define WAVELETTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with wavelet base functions with boundaries, pentagon cut
     */
    class WaveletBoundaryGrid : public Grid {
      protected:
        WaveletBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with wavelet base functions with boundaries, pentagon cut
         *
         * @param dim the dimension of the grid
         */
        WaveletBoundaryGrid(size_t dim);

        /**
         * Destructor
         */
        virtual ~WaveletBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

    };

  }
}

#endif /* WAVELETTRUNCATEDBOUNDARYGRID_HPP */
