/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
#ifndef PREWAVELETGRID_HPP
#define PREWAVELETGRID_HPP

#include "base/grid/Grid.hpp"

namespace sg {
  namespace base {

    /**
     * grid with prewavelet base functions
     */
    class PrewaveletGrid : public Grid {
      protected:
        PrewaveletGrid(std::istream& istr);
        GridStorage* shadowStorage;

      public:
        PrewaveletGrid(size_t dim);
        virtual ~PrewaveletGrid();

        virtual const char* getType();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

        /**
         * gets a pointer to the GridStorage object
         *
         * @return pointer to the GridStorage object
         */
        GridStorage* getShadowStorage();

    };

  }
}

#endif /* PREWAVELETGRID_HPP */
