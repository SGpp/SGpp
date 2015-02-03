// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CLENSHAWCURTISTABLE_HPP
#define CLENSHAWCURTISTABLE_HPP

#include <vector>
#include <cmath>
#include <cstddef>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Lookup table for 1D Clenshaw-Curtis points.
     * This class precomputes the first \c maxLevel levels of a 1D Clenshaw-Curtis
     * grid to increase performance of Clenshaw-Curtis grids.
     */
    template<class LT, class IT>
    class ClenshawCurtisTable {
      public:
        typedef LT level_type;
        typedef IT index_type;

        /// default number of intervals
        static const LT DEFAULT_MAX_LEVEL = 16;

        /**
         * Constructor creating the lookup table.
         *
         * @param maxLevel    level up to which grid points should be pre-computed
         */
        ClenshawCurtisTable(LT maxLevel = DEFAULT_MAX_LEVEL);

        /**
         * @param l       level of the grid point
         * @param i       index of the grid point (can be even)
         */
        inline double getPoint(LT l, IT i) const {
          if (l <= maxLevel) {
            return table[(1 << l) + l + i - 1];
          } else {
            const double h = 1.0 / static_cast<double>(1 << l);
            return (cos(M_PI * (1.0 - static_cast<double>(i) * h)) + 1.0) / 2.0;
          }
        }

      protected:
        /// lookup table
        std::vector<double> table;
        /// maximal level
        LT maxLevel;
    };

    /// typedef for standard level/index types
    typedef ClenshawCurtisTable<unsigned int, unsigned int> SClenshawCurtisTable;

    extern SClenshawCurtisTable clenshawCurtisTable;

  }
}

#endif /* CLENSHAWCURTISTABLE_HPP */
