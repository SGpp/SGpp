// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace base {

/**
 * Polynomial basis functions.
 *
 * @version $HEAD$
 */
template <class LT, class IT>
class PolyClenshawCurtisBasis : public Basis<LT, IT> {
 public:
  /**
   * Constructor
   *
   * @param degree the polynom's max. degree
   */
  explicit PolyClenshawCurtisBasis(size_t degree)
      : degree(degree),
        clenshawCurtisTable(ClenshawCurtisTable::getInstance()),
        idxtable(4),
        quadRule(GaussLegendreQuadRule1D::getInstance()) {
    if (degree < 2) {
      throw factory_exception("PolyClensahwCurtisBasis: degree < 2");
    }

    if (degree > 20) {
      throw factory_exception("PolyClensahwCurtisBasis: degree > 20 is not supported");
    }

    idxtable[0] = 1;
    idxtable[1] = 2;
    idxtable[2] = -2;
    idxtable[3] = -1;
  }

  /**
   * Destructor
   */
  ~PolyClenshawCurtisBasis() override {}

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

  size_t getDegree() const override { return degree; }

  double evalDx(LT level, IT index, double x) {
    // uses the logarithmic derivative method from the second answer
    // http://math.stackexchange.com/questions/809927/first-derivative-of-lagrange-polynomial

    size_t deg = std::min<size_t>(degree, level + 1);
    double result = eval(level, index, x);
    if (result == 0.0) return 0.0;

    double sum = 0.0;
    // see eval-function for explanation of traversal code
    size_t root = index;
    size_t id = root;
    root++;
    HashGridPoint gp(1);
    gp.setAsHierarchicalGridPoint(0, level, static_cast<IT>(root));

    double xroot = clenshawCurtisTable.getPoint(gp.getLevel(0), gp.getIndex(0));
    sum += 1 / (x - xroot);
    root -= 2;
    for (size_t j = 2; j < static_cast<size_t>(1 << deg); j *= 2) {
      gp.setAsHierarchicalGridPoint(0, level, static_cast<IT>(root));
      xroot = clenshawCurtisTable.getPoint(gp.getLevel(0), gp.getIndex(0));
      sum += 1 / (x - xroot);
      root += idxtable[id & 3] * j;
      id >>= 1;
    }
    return result * sum;
  }

  double eval(LT level, IT index, double p) override {
    // check if x value is in the unit interval
    if ((p <= 0) || (p >= 1)) {
      return 0.0f;
    }

    // load boundaries of support
    double xleft = clenshawCurtisTable.getPoint(level, index - 1);
    double xright = clenshawCurtisTable.getPoint(level, index + 1);

    // check if p is out of bounds
    if ((p <= xleft) || (p >= xright)) {
      return 0.0f;
    } else {
      return evalBasis(level, index, p);
    }
  }

  /**
   * Evaluate a basis function.
   * Has a dependence on the absolute position of grid point and support.
   *
   * We compute the roots in units of h = grid spacing at level l = 2 ** -l.
   *
   * Due to limited polynomial degree, we compute the roots of the Lagrange
   * polynomial bottom up.
   */
  double evalBasis(LT level, IT index, double p) {
    // degree of polynomial, limited with level of grid point
    size_t deg = std::min<size_t>(degree, level + 1);
    // start with the current grid point
    IT root = index;
    // copy of index: used to identify the path in the binary tree of grid
    // points. The binary representation of the index contains the information
    // in which direction the grid point has been added w.r.t.
    // the parent node.
    // 0011 -> left, left, right (last one ignored)
    IT id = root;
    // position where the polynomial is one -> position of grid point:
    // (level, index)
    double xbase = clenshawCurtisTable.getPoint(level, index);
    double eval = 1;
    // as first parent we choose the right one. In units of h, it is 1 distance
    // away from the current one.
    root++;
    HashGridPoint gp(1);
    gp.setAsHierarchicalGridPoint(0, level, root);

    // add it to the Lagrange polynomial and normalize it
    double xroot = clenshawCurtisTable.getPoint(gp.getLevel(0), gp.getIndex(0));
    eval *= (p - xroot) / (xbase - xroot);
    // go to the next left neighbor that must exist due to
    // minimum degree of 2 of
    // the polynomial. the reference point is now the last one
    // stored in root, which
    // is the right neighbor of p. So here we need to go 2 units
    // of h to the left.
    root -= 2;

    // p - 1 runs in this loop: so in total the polynomial has a
    // degree of p taking
    // into account that the first root has been added already
    for (unsigned int j = 2; j < static_cast<unsigned int>(1 << deg); j *= 2) {
      // add the next root to the polynomial
      gp.setAsHierarchicalGridPoint(0, level, root);
      xroot = clenshawCurtisTable.getPoint(gp.getLevel(0), gp.getIndex(0));
      eval *= (p - xroot) / (xbase - xroot);
      // take last two indices (id & 3):
      // this gives you the information where to
      // look for the next root. The result needs to be scaled
      // with j due to the fact
      // that we calculate the roots in units of h.
      // We go bottom up, therefore the
      // distance to the next root changes by a factor of +- {1, 2}
      // depending of the history
      // of the grid point. We just consider the history in the
      // last two indices, so the
      // relation between the current grid point to its first predecessor.
      // The scaling by j is due to fact that the distance needs to be
      // computed in units of h.
      root += idxtable[id & 3] * j;
      // remove the last index that means that we go one step up in the tree.
      // Then do the same again until we reach the maximum polynomial degree.
      id >>= 1;
    }

    return eval;
  }

  double getIntegral(LT level, IT index) {
    // load boundaries of support
    double xleft = clenshawCurtisTable.getPoint(level, index - 1);
    double xright = clenshawCurtisTable.getPoint(level, index + 1);

    // grid spacing
    double h = xright - xleft;

    // --------------------------------
    // Gauss-Legendre quadrature
    // --------------------------------
    size_t deg = std::min<size_t>(degree, level + 1);
    size_t n_roots = ((deg + 1) >> 1) + 1;  // ceil((deg + 1) / 2) - 1
    base::DataVector roots(n_roots);
    base::DataVector weights(n_roots);
    // getting legendre gauss points and weights in [-1, 1]
    quadRule.getLevelPointsAndWeights(n_roots, roots, weights);

    double sum = 0.0f;
    double x = 0.0f;

    for (size_t i = 0; i < n_roots; i++) {
      // scale the roots to the support of the basis:
      // [-1, 1] -> [0, 1] -> [a, b]
      x = (1 + roots[i]) / 2. * h + xleft;
      // evaluate the polynom and weight it
      sum += weights[i] * eval(level, index, x);
    }

    // scale the result with the width of the support
    return h * sum / 2.;
  }

 protected:
  /// the polynom's max degree
  size_t degree;
  /// reference to the Clenshaw-Curtis cache table
  ClenshawCurtisTable& clenshawCurtisTable;
  // compute values for roots
  std::vector<int> idxtable;

 private:
  /// gauss legendre quadrature rule to compute the integral of the bases
  base::GaussLegendreQuadRule1D& quadRule;
};

// default type-def (unsigned int for level and index)
typedef PolyClenshawCurtisBasis<unsigned int, unsigned int> SPolyClenshawCurtisBase;

}  // namespace base
}  // namespace sgpp
