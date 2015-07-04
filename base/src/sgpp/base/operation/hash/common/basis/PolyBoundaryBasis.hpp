// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POLYBOUNDARY_BASE_HPP
#define POLYBOUNDARY_BASE_HPP

#include <cmath>
#include <vector>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Polynomial basis functions with boundaries.
     *
     * @version $HEAD$
     */
    template<class LT, class IT>
    class PolyBoundaryBasis: public Basis<LT, IT> {
      public:
        /**
         * Constructor
         *
         * @param degree the polynom's max. degree
         */
        PolyBoundaryBasis(size_t degree) :
          degree(degree) {
          if (degree < 2) {
            throw factory_exception("PolyBoundaryBasis: degree < 2");
          }

          if (degree > 20) {
            throw factory_exception(
              "PolyBoundaryBasis: degree > 20 is not supported");
          }

          idxtable = new int[4];
          idxtable[0] = 1;
          idxtable[1] = 2;
          idxtable[2] = -2;
          idxtable[3] = -1;
        }

        /**
         * Destructor
         */
        ~PolyBoundaryBasis() {
          delete idxtable;
        }

        float_t evalHierToTop(LT level, IT index, DataVector& coeffs, float_t pos) {
          float_t result = 0.0;

          if (level > 0) {
            // add left and right boundary point
            result += coeffs[0] * eval(0, 0, pos) + coeffs[1] * eval(0, 1, pos);

            // just evaluate the hierarchical ancestors -> so start with the
            // parent node
            level--;
            index >>= 1;
            index |= 1;
          }

          for (; level >= 1; level--) {
            // the coefficients are shifted by one due to the two boundary points
            // ergo: coeffs = [<left boundary, <right boundary>, 1, 2, ...]
            result += coeffs[level + 1] * eval(level, index, pos);
            index >>= 1; // == ((index - 1) / 2)
            index |= 1; // == (index % 2 == 0) ? (index + 1) : index;
          }

          return result;
        }

        size_t getDegree() {
          return degree;
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
        float_t eval(LT level, IT index, float_t p) {
          if (level == 0) {
            if (index == 0) {
              return 1 - p;
            } else {
              // index == 1
              return p;
            }
          }

          // degree of polynomial, limited with level of grid point
          size_t deg = std::min<size_t>(degree, level + 1);
          // get the position in units of h of the current maximum level
          p *= (1 << level);
          // start with the current grid point
          size_t root = index;
          // copy of index: used to identify the path in the binary tree of grid
          // points. The binary representation of the index contains the information
          // in which direction the grid point has been added w.r.t. to the parent node.
          // 0011 -> left, left, right (last one ignored)
          size_t id = root;
          // position where the polynomial is one -> position of grid point: (level, index)
          float_t base = static_cast<float_t>(root);
          float_t eval = 1;
          // as first parent we choose the right one. In units of h, it is 1 distance
          // away from the current one.
          root++;
          // add it to the Lagrange polynomial and normalize it
          eval *= (p - root) / (base - root);
          // go to the next left neighbor that must exist due to minimum degree of 2 of
          // the polynomial. the reference point is now the last one stored in root, which
          // is the right neighbor of p. So here we need to go 2 units of h to the left.
          root -= 2;

          // p - 1 runs in this loop: so in total the polynomial has a degree of p taking
          // into account that the first root has been added already
          for (size_t j = 2; j < (1 << deg); j *= 2) {
            // add the next root to the polynomial
            eval *= (p - root) / (base - root);
            // take last two indices (id & 3): this gives you the information where to
            // look for the next root. The result needs to be scaled with j due to the fact
            // that we calculate the roots in units of h. We go bottom up, therefore the
            // distance to the next root changes by a factor of +- {1, 2} depending of the history
            // of the grid point. We just consider the history in the last two indices, so the
            // relation between the current grid point to its first predecessor.
            // The scaling by j is due to fact that the distance needs to be
            // computed in units of h.
            root += idxtable[id & 3] * j;
            // remove the last index that means that we go one step up in the tree. Then do
            // the same again until we reach the maximum polynomial degree.
            id >>= 1;
          }

          return eval;
        }

        float_t evalSave(LT level, IT index, float_t p) {
          // spacing on current level
          float_t h = 1.0f / (1 << level);

          // check if p is out of bounds
          if (p <= h * (index - 1) || p >= h * (index + 1))
            return 0.0f;
          else
            return eval(level, index, p);
        }

        float_t eval(LT level, IT index, float_t p, float_t offset, float_t width) {
          // for bounding box evaluation
          // scale p in [offset, offset + width] linearly to [0, 1] and do simple
          // evaluation
          return eval(level, index, (p - offset) / width);
        }

        float_t getIntegral(LT level, IT index) {
          if (level == 0) {
            return 0.5;
          }

          // grid spacing
          float_t h = 1.0f / (1 << level);

          // --------------------------------
          // Gauss-Legendre quadrature
          // --------------------------------
          size_t deg = std::min<size_t>(degree, level + 1);
          size_t n_roots = ((deg + 1) >> 1) + 1;        // ceil((deg + 1) / 2) - 1
          base::DataVector roots(n_roots);
          base::DataVector weights(n_roots);
          // getting legendre gauss points and weights in [-1, 1]
          quadRule.getLevelPointsAndWeights(n_roots, roots, weights);

          float_t sum = 0.0f;
          float_t x = 0.0f;

          for (int i = 0; i < n_roots; i++) {
            // scale the roots to the support of the basis:
            // [-1, 1] -> [0, 1] -> [a, b]
            x = h * (roots[i] + index);
            // evaluate the polynom and weight it
            sum += weights[i] * eval(level, index, x);
          }

          // scale the result with the width of the support
          return h * sum;
        }

      protected:
        /// the polynom's max degree
        size_t degree;
        /// compute values for roots
        int* idxtable;
      private:
        /// gauss legendre quadrature rule to compute the integral of the bases
        base::GaussLegendreQuadRule1D quadRule;
    };

    // default type-def (unsigned int for level and index)
    typedef PolyBoundaryBasis<unsigned int, unsigned int> SPolyBoundaryBase;

  }
}

#endif /* POLYBOUNDARY_BASE_HPP */
