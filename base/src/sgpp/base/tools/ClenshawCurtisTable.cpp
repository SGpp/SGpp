// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

namespace SGPP {
  namespace base {

    template<class LT, class IT>
    ClenshawCurtisTable<LT, IT>::ClenshawCurtisTable(LT maxLevel)
      : table((1 << (maxLevel + 1)) + maxLevel),
        maxLevel(maxLevel) {
      size_t k = 0;
      IT hInv = 1;

      for (LT l = 0; l <= maxLevel; l++) {
        const float_t h = 1.0 / static_cast<float_t>(hInv);

        for (IT i = 0; i <= hInv; i++) {
          table[k] = (cos(M_PI * (1.0 - static_cast<float_t>(i) * h)) + 1.0) / 2.0;
          k++;
        }

        hInv *= 2;
      }
    }

    SClenshawCurtisTable clenshawCurtisTable;

  }
}
