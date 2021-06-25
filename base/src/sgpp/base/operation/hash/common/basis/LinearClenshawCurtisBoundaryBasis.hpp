// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * Linear basis on Clenshaw-Curtis grids.
 */
template <class LT, class IT>
class LinearClenshawCurtisBoundaryBasis : public Basis<LT, IT> {
 public:
  LinearClenshawCurtisBoundaryBasis() {}

  /**
   * Destructor.
   */
  ~LinearClenshawCurtisBoundaryBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of Clenshaw-Curtis linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 0) {
      // first level
      if (i == 0) {
        return 1.0 - x;
      } else {
        return x;
      }
    } else {
      return basis.eval(l, i, x);
    }
  }

  double eval(LT level, IT index, double p, double offset, double width) {
    // for bounding box evaluation
    // scale p in [offset, offset + width] linearly to [0, 1] and do simple
    // evaluation
    return eval(level, index, (p - offset) / width);
  }
  double evalDx(LT level, IT index, double x) override {
    std::cerr << "LinearClenshawCurtisBoundaryBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  double getIntegral(LT level, IT index) override {
    // boundary points
    if (level == 0) {
      return 0.5;
    } else {
      return basis.getIntegral(level, index);
    }
  }

  inline size_t getDegree() const override { return 1; }

 protected:
  /// linear clenshaw curtis basis
  SLinearClenshawCurtisBase basis;
};

// default type-def (unsigned int for level and index)
typedef LinearClenshawCurtisBoundaryBasis<unsigned int, unsigned int>
    SLinearClenshawCurtisBoundaryBase;

}  // namespace base
}  // namespace sgpp
