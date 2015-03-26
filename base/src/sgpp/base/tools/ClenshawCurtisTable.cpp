// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

namespace SGPP {
  namespace base {

    ClenshawCurtisTable::ClenshawCurtisTable(level_type maxLevel)
      : table((1 << (maxLevel + 1)) + maxLevel),
        maxLevel(maxLevel) {
      size_t k = 0;
      index_type hInv = 1;

      for (level_type l = 0; l <= maxLevel; l++) {
        const float_t h = 1.0 / static_cast<float_t>(hInv);

        for (index_type i = 0; i <= hInv; i++) {
          table[k] = calculatePoint(h, i);
          k++;
        }

        hInv *= 2;
      }
    }

    ClenshawCurtisTable clenshawCurtisTable;

  }
}
