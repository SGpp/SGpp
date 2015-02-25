// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PREWAVELETGRID_HPP
#define PREWAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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

        virtual const SBasis& getBasis();

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