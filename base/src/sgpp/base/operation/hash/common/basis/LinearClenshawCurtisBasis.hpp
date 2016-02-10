// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_CLENSHAW_CURTIS_BASE_HPP
#define LINEAR_CLENSHAW_CURTIS_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace SGPP {
namespace base {

/**
 * Linear basis on Clenshaw-Curtis grids.
 */
template<class LT, class IT>
class LinearClenshawCurtisBasis: public Basis<LT, IT> {
 public:
  LinearClenshawCurtisBasis() :
    clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {
  }

  /**
   * Destructor.
   */
  ~LinearClenshawCurtisBasis() override {
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of Clenshaw-Curtis linear basis function
   */
  inline float_t eval(LT l, IT i, float_t x) override {
    if (l == 0) {
      // first level
      if (i == 0) {
        return 1.0 - x;
      } else {
        return x;
      }
    } else {
      // endpoints of support
      const float_t x0 = clenshawCurtisTable.getPoint(l, i - 1);
      const float_t x2 = clenshawCurtisTable.getPoint(l, i + 1);

      if ((x <= x0) || (x >= x2)) {
        // point out of support
        return 0.0;
      }

      // peak of basis function
      const float_t x1 = clenshawCurtisTable.getPoint(l, i);

      // linear interpolation between (x0, x1, x2), (0, 1, 0)
      if (x < x1) {
        return 1.0 - (x1 - x) / (x1 - x0);
      } else {
        return (x2 - x) / (x2 - x1);
      }
    }
  }

 protected:
  /// reference to the Clenshaw-Curtis cache table
  ClenshawCurtisTable& clenshawCurtisTable;
};

// default type-def (unsigned int for level and index)
typedef LinearClenshawCurtisBasis<unsigned int, unsigned int>
SLinearClenshawCurtisBase;

}  // namespace base
}  // namespace SGPP

#endif /* LINEAR_CLENSHAW_CURTIS_BASE_HPP */
