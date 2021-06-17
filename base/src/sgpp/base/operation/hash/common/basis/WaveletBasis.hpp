// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELET_BASE_HPP
#define WAVELET_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * Wavelet basis on Noboundary grids.
 */
template <class LT, class IT>
class WaveletBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~WaveletBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of wavelet basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);
    const double t = x * hInv - static_cast<double>(i);

    if ((t >= 2.0) || (t <= -2.0)) {
      // out of support (cut-off)
      return 0.0;
    }

    const double t2 = t * t;
    return (1.0 - t2) * std::exp(-t2);
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of wavelet basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);
    const double t = x * hInv - static_cast<double>(i);

    if ((t >= 2.0) || (t <= -2.0)) {
      // out of support (cut-off)
      return 0.0;
    }

    const double t2 = t * t;
    return 2.0 * t * (t2 - 2.0) * std::exp(-t2) * hInv;
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of wavelet basis function
   */
  inline double evalDxDx(LT l, IT i, double x) {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);
    const double t = x * hInv - static_cast<double>(i);

    if ((t >= 2.0) || (t <= -2.0)) {
      // out of support (cut-off)
      return 0.0;
    }

    const double t2 = t * t;
    return -2.0 * (2.0 * t2 * t2 - 7.0 * t2 + 2.0) * std::exp(-t2) * hInv * hInv;
  }

  inline double getIntegral(LT level, IT index) override { return -1.0; }

  inline size_t getDegree() const override { return 0; }
};

// default type-def (unsigned int for level and index)
typedef WaveletBasis<unsigned int, unsigned int> SWaveletBase;

}  // namespace base
}  // namespace sgpp

#endif /* WAVELET_BASE_HPP */
