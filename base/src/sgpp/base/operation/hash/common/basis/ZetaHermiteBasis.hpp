// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * Zeta basis function for Hermite Interpolation
    f(0)=0
    f'(0)=1
 */
template <class LT, class IT>
class ZetaHermiteBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~ZetaHermiteBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of psi basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);

    double x_uniform = x * hInv - static_cast<double>(i);

    if (x_uniform > 1 || x_uniform < -1)
      return 0;

    else if (x_uniform  < 0) {
      return pow(x_uniform, 3) + 2 * pow(x_uniform, 2) + x_uniform;
      
      
    } else {
      return pow(x_uniform, 3) - 2 * pow(x_uniform, 2) + x_uniform;
      
    }
  }
};

}  // namespace base
}  // namespace sgpp
