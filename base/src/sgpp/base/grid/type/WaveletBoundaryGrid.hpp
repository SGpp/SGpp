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
         * @param boundaryLevel 1 + how much levels the boundary is coarser than
         *                      the main axes, 0 means one level finer,
         *                      1 means same level,
         *                      2 means one level coarser, etc.
         */
        WaveletBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

        /**
         * Destructor
         */
        virtual ~WaveletBoundaryGrid() override;

        virtual SGPP::base::GridType getType() override;

        virtual const SBasis& getBasis() override;

        virtual GridGenerator* createGridGenerator() override;

        static Grid* unserialize(std::istream& istr);

      protected:
        /// 1 + how much levels the boundary is coarser than the main axes
        level_t boundaryLevel;
    };

  }
}

#endif /* WAVELETTRUNCATEDBOUNDARYGRID_HPP */
