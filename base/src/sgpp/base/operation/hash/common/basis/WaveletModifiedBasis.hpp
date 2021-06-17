// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELET_MODIFIED_BASE_HPP
#define WAVELET_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * Modified wavelet basis on Noboundary grids.
 */
template <class LT, class IT>
class WaveletModifiedBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~WaveletModifiedBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of modified wavelet basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 1) {
      // first level
      return 1.0;
    }

    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    const double h = 1.0 / hInvDbl;

    if ((i == 1) && (x < 1.560231504260063 * h)) {
      // left modified basis function
      return 0.501309319347014 + 1.38033323862282 * (0.560231504260063 - x * hInvDbl + 1.0);
    } else if ((i == hInv - 1) && (x > 1.0 - 1.560231504260063 * h)) {
      // right modified basis function
      return 0.501309319347014 +
             1.38033323862282 * (0.560231504260063 + x * hInvDbl - static_cast<double>(i));
    } else {
      // interior basis function
      // (or modified function, but x is in the non-modified part)
      const double t = x * hInvDbl - static_cast<double>(i);

      if ((t > 2.0) || (t < -2.0)) {
        // out of support (cut-off)
        return 0.0;
      }

      const double t2 = t * t;
      return (1.0 - t2) * std::exp(-t2);
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of modified wavelet basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    if (l == 1) {
      // first level
      return 0.0;
    }

    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    const double h = 1.0 / hInvDbl;

    if ((i == 1) && (x < 1.560231504260063 * h)) {
      // left modified basis function
      return -1.38033323862282 * hInvDbl;
    } else if ((i == hInv - 1) && (x > 1.0 - 1.560231504260063 * h)) {
      // right modified basis function
      return 1.38033323862282 * hInvDbl;
    } else {
      // interior basis function
      // (or modified function, but x is in the non-modified part)
      const double t = x * hInvDbl - static_cast<double>(i);

      if ((t > 2.0) || (t < -2.0)) {
        // out of support (cut-off)
        return 0.0;
      }

      const double t2 = t * t;
      return 2.0 * t * (t2 - 2.0) * std::exp(-t2) * hInvDbl;
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of modified wavelet basis function
   */
  inline double evalDxDx(LT l, IT i, double x) {
    if (l == 1) {
      // first level
      return 0.0;
    }

    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    const double h = 1.0 / hInvDbl;

    if ((i == 1) && (x < 1.560231504260063 * h)) {
      // left modified basis function
      return 0.0;
    } else if ((i == hInv - 1) && (x > 1.0 - 1.560231504260063 * h)) {
      // right modified basis function
      return 0.0;
    } else {
      // interior basis function
      // (or modified function, but x is in the non-modified part)
      const double t = x * hInvDbl - static_cast<double>(i);

      if ((t > 2.0) || (t < -2.0)) {
        // out of support (cut-off)
        return 0.0;
      }

      const double t2 = t * t;
      const double t4 = t2 * t2;
      return -2.0 * (2.0 * t4 - 7.0 * t2 + 2.0) * std::exp(-t2) * hInvDbl * hInvDbl;
    }
  }

  inline double getIntegral(LT level, IT index) override { return -1.0; }

  inline size_t getDegree() const override { return 0; }
};

// default type-def (unsigned int for level and index)
typedef WaveletModifiedBasis<unsigned int, unsigned int> SWaveletModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* WAVELET_MODIFIED_BASE_HPP */
