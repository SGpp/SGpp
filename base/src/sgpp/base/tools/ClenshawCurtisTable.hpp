// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CLENSHAWCURTISTABLE_HPP
#define CLENSHAWCURTISTABLE_HPP

#include <cmath>
#include <cstddef>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>

namespace SGPP {
  namespace base {

    /**
     * Lookup table for 1D Clenshaw-Curtis points.
     * This class precomputes the first \c maxLevel levels of a 1D Clenshaw-Curtis
     * grid to increase performance of Clenshaw-Curtis grids.
     */
    class ClenshawCurtisTable {
      public:
        typedef HashGridIndex::level_type level_type;
        typedef HashGridIndex::index_type index_type;

        /// default number of intervals
        static const level_type DEFAULT_MAX_LEVEL = 16;

        /**
         * Constructor creating the lookup table.
         *
         * @param maxLevel    level up to which grid points should be pre-computed
         */
        ClenshawCurtisTable(level_type maxLevel = DEFAULT_MAX_LEVEL);

        /**
         * @param l       level of the grid point
         * @param i       index of the grid point (can be even)
         */
        inline float_t getPoint(level_type l, index_type i) const {
          if (l <= maxLevel) {
            return table.get((1 << l) + l + i - 1);
          } else {
            const float_t h = 1.0 / static_cast<float_t>(1 << l);
            return (cos(M_PI * (1.0 - static_cast<float_t>(i) * h)) + 1.0) / 2.0;
          }
        }

      protected:
        /// lookup table
        DataVector table;
        /// maximal level
        level_type maxLevel;
    };

    extern ClenshawCurtisTable clenshawCurtisTable;

  }
}

#endif /* CLENSHAWCURTISTABLE_HPP */
