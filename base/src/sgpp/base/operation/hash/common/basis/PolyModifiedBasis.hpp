// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>

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
class PolyModifiedBasis : public Basis<LT, IT> {
 protected:
  /// poly basis
  SPolyBase polyBasis;

 public:
  /**
   * Constructor
   *
   * @param degree the polynom's max. degree
   */
  explicit PolyModifiedBasis(size_t degree) : polyBasis(degree) {}

  /**
   * Destructor
   */
  ~PolyModifiedBasis() override {}

  double eval(LT level, IT index, double p) override {
    // spacing on current level
    double h = 1.0 / static_cast<double>(1 << level);

    // check if p is out of bounds
    if ((p < h * static_cast<double>(index - 1)) || (p > h * static_cast<double>(index + 1))) {
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

  double evalDx(LT level, IT index, double x) override {
    double hInvDbl = static_cast<double>(1 << level);
    const IT hInv = static_cast<IT>(1) << level;
    if (level == 1) {
      // first level
      return 0.0;
    } else if (index == 1) {
      return ((x <= 2.0 / hInvDbl) ? -hInvDbl : 0.0);
    } else if (index == hInv - 1) {
      return ((x >= 1.0 - 2.0 / hInvDbl) ? hInvDbl : 0.0);
    } else {
      // interior basis function
      return polyBasis.evalDx(level, index, x);
    }
  }

  double getIntegral(LT level, IT index) override {
    const IT hInv = static_cast<IT>(1) << level;

    if (level == 1) {
      // first level
      return 1.0;
    } else if ((index == 1) || (index == hInv - 1)) {
      // left and right modified basis functions
      return 2. / static_cast<double>(hInv);
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
      //        index = ((index - 1) / 2);
      //        index = (index % 2 == 0) ? (index + 1) : index;
    }

    return result;
  }

  size_t getDegree() const override { return polyBasis.getDegree(); }

 private:
  /**
   * Evaluate a basis function.
   * Has a dependence on the absolute position of grid point and support.
   */
  double evalBasis(LT level, IT index, double p) {
    const IT hInv = static_cast<IT>(1) << level;
    const double hInvDbl = static_cast<double>(hInv);

    if (level == 1) {
      // first level
      return 1.0;
    } else if (index == 1) {
      // left modified basis function
      return ((p <= 2.0 / hInvDbl) ? (2.0 - hInvDbl * p) : 0.0);
    } else if (index == hInv - 1) {
      // right modified basis function
      return ((p >= 1.0 - 2.0 / hInvDbl) ? (hInvDbl * p - static_cast<double>(index) + 1.0) : 0.0);
    } else {
      // interior basis function
      return polyBasis.eval(level, index, p);
    }
  }
};

// default type-def (unsigned int for level and index)
typedef PolyModifiedBasis<unsigned int, unsigned int> SPolyModifiedBase;

}  // namespace base
}  // namespace sgpp
