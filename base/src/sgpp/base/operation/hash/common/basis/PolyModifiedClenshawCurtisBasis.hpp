// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <cmath>
#include <vector>

namespace sgpp {
namespace base {

/**
 * Modified polynomial base functions.
 * Special polynomial functions to cover values unequal 0 at the border. Implemented as seen in AWR
 * 2 paper
 * by Prof. Bungartz
 * (http://www5.in.tum.de/wiki/index.php/Algorithmen_des_Wissenschaftlichen_Rechnens_II_-_Winter_08)
 */
template <class LT, class IT>
class PolyModifiedClenshawCurtisBasis : public Basis<LT, IT> {
 protected:
  /// poly basis
  SPolyClenshawCurtisBase polyBasis;
  /// reference to the Clenshaw-Curtis cache table
  ClenshawCurtisTable& clenshawCurtisTable;

 public:
  /**
   * Constructor
   *
   * @param degree the polynom's max. degree
   */
  explicit PolyModifiedClenshawCurtisBasis(size_t degree)
      : polyBasis(degree), clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {}

  /**
   * Destructor
   */
  ~PolyModifiedClenshawCurtisBasis() override {}

  double eval(LT level, IT index, double p) override {
    // load boundaries of support
    double xleft = clenshawCurtisTable.getPoint(level, index - 1);
    double xright = clenshawCurtisTable.getPoint(level, index + 1);

    // check if p is out of bounds
    if ((p < xleft) || (p > xright)) {
      return 0.0;
    } else {
      return evalBasis(level, index, p);
    }
  }

  double eval(LT level, IT index, double p, double offset, double width) {
    // for bounding box evaluation
    // scale p in [offset, offset + width] linearly to [0, 1] and do simple
    // evaluation
    return eval(level, index, (p - offset) / width);
  }

  double getIntegral(LT level, IT index) override {
    const IT hInv = static_cast<IT>(1) << level;

    if (level == 1) {
      // first level
      return 1.0;
    } else if ((index == 1) || (index == hInv - 1)) {
      // left and right modified basis functions
      // load right boundary of support
      double xr = clenshawCurtisTable.getPoint(level, 2);
      double xc = clenshawCurtisTable.getPoint(level, 1);
      return 0.5 * (xc / (xr - xc) + 1.0) * xr;
    } else {
      // interior basis function
      return polyBasis.getIntegral(level, index);
    }
  }

  /**
   * Evaluates all the hierarchical ancestors of the node defined by level
   * and index. NOTE: It does not evaluate the current node itself.
   *
   * @param level
   * @param index
   * @param coeffs
   * @param pos
   * @return
   */
  double evalHierToTop(LT level, IT index, DataVector& coeffs, double pos) {
    double result = 0.0;

    // just evaluate the hierarchical ancestors -> so start with the
    // parent node
    level--;
    index >>= 1;
    index |= 1;

    for (; level >= 1; level--) {
      result += coeffs[level] * eval(level, index, pos);
      index >>= 1;
      index |= 1;
    }

    return result;
  }

  size_t getDegree() const override { return polyBasis.getDegree(); }

  double evalDx(LT level, IT index, double x) override {
    const IT hInv = static_cast<IT>(1) << level;
    if (level == 1) {
      // first level
      return 0.0;
    } else if (index == 1) {
      // left modified basis function
      double xr = clenshawCurtisTable.getPoint(level, index + 1);
      double xc = clenshawCurtisTable.getPoint(level, index);
      return ((x < xr) ? (-1.0 / (xr - xc)) : 0.0);
    } else if (index == hInv - 1) {
      // right modified basis function
      double xl = clenshawCurtisTable.getPoint(level, index - 1);
      double xc = clenshawCurtisTable.getPoint(level, index);
      return ((x > xl) ? (1.0 / (xc - xl)) : 0.0);
    } else {
      // interior basis function
      return polyBasis.eval(level, index, x);
    }
  }

 private:
  /**
   * Evaluate a basis function.
   * Has a dependence on the absolute position of grid point and support.
   */
  double evalBasis(LT level, IT index, double p) {
    const IT hInv = static_cast<IT>(1) << level;
    if (level == 1) {
      // first level
      return 1.0;
    } else if (index == 1) {
      // left modified basis function
      double xr = clenshawCurtisTable.getPoint(level, index + 1);
      double xc = clenshawCurtisTable.getPoint(level, index);
      return ((p < xr) ? (-1.0 / (xr - xc) * (p - xc) + 1.0) : 0.0);
    } else if (index == hInv - 1) {
      // right modified basis function
      double xl = clenshawCurtisTable.getPoint(level, index - 1);
      double xc = clenshawCurtisTable.getPoint(level, index);
      return ((p > xl) ? (1.0 / (xc - xl) * (p - xl)) : 0.0);
    } else {
      // interior basis function
      return polyBasis.eval(level, index, p);
    }
  }
};

// default type-def (unsigned int for level and index)
typedef PolyModifiedClenshawCurtisBasis<unsigned int, unsigned int> SPolyModifiedClenshawCurtisBase;

}  // namespace base
}  // namespace sgpp
