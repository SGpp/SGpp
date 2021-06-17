// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * Linear basis on Clenshaw-Curtis grids.
 */
template <class LT, class IT>
class LinearClenshawCurtisBasis : public Basis<LT, IT> {
 public:
  LinearClenshawCurtisBasis() : clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {}

  /**
   * Destructor.
   */
  ~LinearClenshawCurtisBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of Clenshaw-Curtis linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    // endpoints of support
    const double x0 = clenshawCurtisTable.getPoint(l, i - 1);
    const double x2 = clenshawCurtisTable.getPoint(l, i + 1);
    // peak of basis function
    const double x1 = clenshawCurtisTable.getPoint(l, i);

    // linear interpolation between (x0, x1, x2), (0, 1, 0)
    if (x < x1) {
      return std::max(0.0, 1.0 - (x1 - x) / (x1 - x0));
    } else {
      return std::max(0.0, (x2 - x) / (x2 - x1));
    }
  }

  double eval(LT level, IT index, double p, double offset, double width) {
    // for bounding box evaluation
    // scale p in [offset, offset + width] linearly to [0, 1] and do simple
    // evaluation
    return eval(level, index, (p - offset) / width);
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "LinearClenshawCurtisBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  double getIntegral(LT level, IT index) override {
    // endpoints of support
    const double x0 = clenshawCurtisTable.getPoint(level, index - 1);
    const double x2 = clenshawCurtisTable.getPoint(level, index + 1);
    return (x2 - x0) / 2.;
  }

  inline size_t getDegree() const override { return 1; }

 protected:
  /// reference to the Clenshaw-Curtis cache table
  ClenshawCurtisTable& clenshawCurtisTable;
};

// default type-def (unsigned int for level and index)
typedef LinearClenshawCurtisBasis<unsigned int, unsigned int> SLinearClenshawCurtisBase;

}  // namespace base
}  // namespace sgpp
