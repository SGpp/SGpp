// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LAGRANGE_SPLINE_BASE_DERIV3_HPP
#define LAGRANGE_SPLINE_BASE_DERIV3_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>
#include <algorithm>

namespace sgpp {
namespace base {

/**
 * Lagrange spline basis (3rd derivative).
 */
template <class LT, class IT>
class LagrangeSplineBasisDeriv3: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  LagrangeSplineBasisDeriv3() : degree(0) {
  }

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit LagrangeSplineBasisDeriv3(size_t degree) : degree(degree) {
    if (degree < 1) {
      this->degree = 1;
    } else if (degree % 2 == 0) {
      this->degree = degree - 1;
    }

    if (this->degree > 7) {
      throw std::runtime_error("Unsupported Lagrange spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~LagrangeSplineBasisDeriv3() override {
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    const double innerDeriv = hInvDbl * hInvDbl * hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (degree) {
      case 1:
        return 0.0;

      case 3:
        if ((t < -3.0) || (t > 3.0)) {
          return 0.0;
        } else if (t < -2.0) {
          return innerDeriv * -2.5000000000000000e-01;
        } else if (t < -1.0) {
          return innerDeriv * 1.7500000000000000e+00;
        } else if (t < 0.0) {
          return innerDeriv * -4.0;
        } else if (t < 1.0) {
          return innerDeriv * 4.0;
        } else if (t < 2.0) {
          return innerDeriv * -1.7500000000000000e+00;
        } else {
          return innerDeriv * 2.5000000000000000e-01;
        }

      case 5:
        if ((t < -5.0) || (t > 5.0)) {
          return 0.0;
        } else if (t < -4.0) {
          t += 5.0;
          double result = 7.5757575757575760e-03;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = -2.3484848484848486e-01;
          result = 1.5151515151515152e-02 + result * t;
          result = 7.5757575757575760e-03 + result * t;
          return innerDeriv * result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = 1.5606060606060606e+00;
          result = -4.5454545454545453e-01 + result * t;
          result = -2.1212121212121213e-01 + result * t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = -4.7424242424242422e+00;
          result = 2.6666666666666665e+00 + result * t;
          result = 8.9393939393939392e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = 8.0;
          result = -6.8181818181818183e+00 + result * t;
          result = -1.1818181818181819e+00 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = -8.0;
          result = 9.1818181818181817e+00 + result * t;
          result *= t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = 4.7424242424242422e+00;
          result = -6.8181818181818183e+00 + result * t;
          result = 1.1818181818181819e+00 + result * t;
          return innerDeriv * result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = -1.5606060606060606e+00;
          result = 2.6666666666666665e+00 + result * t;
          result = -8.9393939393939392e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = 2.3484848484848486e-01;
          result = -4.5454545454545453e-01 + result * t;
          result = 2.1212121212121213e-01 + result * t;
          return innerDeriv * result;
        } else {
          t -= 4.0;
          double result = -7.5757575757575760e-03;
          result = 1.5151515151515152e-02 + result * t;
          result = -7.5757575757575760e-03 + result * t;
          return innerDeriv * result;
        }

      case 7:
        if ((t < -7.0) || (t > 7.0)) {
          return 0.0;
        } else if (t < -6.0) {
          t += 7.0;
          double result = -1.7246136865342163e-05;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -5.0) {
          t += 6.0;
          double result = 2.1902593818984547e-03;
          result = -6.8984547461368653e-05 + result * t;
          result = -1.0347682119205298e-04 + result * t;
          result = -6.8984547461368653e-05 + result * t;
          result = -1.7246136865342163e-05 + result * t;
          return innerDeriv * result;
        } else if (t < -4.0) {
          t += 5.0;
          double result = -3.5389072847682120e-02;
          result = 8.6920529801324496e-03 + result * t;
          result = 1.2831125827814569e-02 + result * t;
          result = 8.2781456953642391e-03 + result * t;
          result = 1.9315673289183224e-03 + result * t;
          return innerDeriv * result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = 2.2951158940397351e-01;
          result = -1.3286423841059603e-01 + result * t;
          result = -1.7342715231788081e-01 + result * t;
          result = -8.1539735099337748e-02 + result * t;
          result = -3.6561810154525387e-03 + result * t;
          return innerDeriv * result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -8.1658733443708609e-01;
          result = 7.8518211920529801e-01 + result * t;
          result = 8.0504966887417218e-01 + result * t;
          result = 9.1059602649006616e-02 + result * t;
          result = -1.6197571743929359e-01 + result * t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 1.8125517384105960e+00;
          result = -2.4811672185430464e+00 + result * t;
          result = -1.7389279801324504e+00 + result * t;
          result = 7.9035596026490063e-01 + result * t;
          result = 7.0272833885209718e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -2.6666666666666665e+00;
          result = 4.7690397350993381e+00 + result * t;
          result = 1.6928807947019868e+00 + result * t;
          result = -2.8807947019867548e+00 + result * t;
          result = -9.1445916114790282e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = 2.6666666666666665e+00;
          result = -5.8976269315673289e+00 + result * t;
          result *= t;
          result = 4.1454194260485648e+00 + result * t;
          result *= t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -1.8125517384105960e+00;
          result = 4.7690397350993381e+00 + result * t;
          result = -1.6928807947019868e+00 + result * t;
          result = -2.8807947019867548e+00 + result * t;
          result = 9.1445916114790282e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = 8.1658733443708609e-01;
          result = -2.4811672185430464e+00 + result * t;
          result = 1.7389279801324504e+00 + result * t;
          result = 7.9035596026490063e-01 + result * t;
          result = -7.0272833885209718e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = -2.2951158940397351e-01;
          result = 7.8518211920529801e-01 + result * t;
          result = -8.0504966887417218e-01 + result * t;
          result = 9.1059602649006616e-02 + result * t;
          result = 1.6197571743929359e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 5.0) {
          t -= 4.0;
          double result = 3.5389072847682120e-02;
          result = -1.3286423841059603e-01 + result * t;
          result = 1.7342715231788081e-01 + result * t;
          result = -8.1539735099337748e-02 + result * t;
          result = 3.6561810154525387e-03 + result * t;
          return innerDeriv * result;
        } else if (t < 6.0) {
          t -= 5.0;
          double result = -2.1902593818984547e-03;
          result = 8.6920529801324496e-03 + result * t;
          result = -1.2831125827814569e-02 + result * t;
          result = 8.2781456953642391e-03 + result * t;
          result = -1.9315673289183224e-03 + result * t;
          return innerDeriv * result;
        } else {
          t -= 6.0;
          double result = 1.7246136865342163e-05;
          result = -6.8984547461368653e-05 + result * t;
          result = 1.0347682119205298e-04 + result * t;
          result = -6.8984547461368653e-05 + result * t;
          result = 1.7246136865342163e-05 + result * t;
          return innerDeriv * result;
        }

      default:
        return 0.0;
    }
  }

  /**
   * @return      Spline degree
   */
  inline size_t getDegree() const {
    return degree;
  }

 protected:
  /// degree of the spline
  size_t degree;
};

// default type-def (unsigned int for level and index)
typedef LagrangeSplineBasisDeriv3<unsigned int, unsigned int> SLagrangeSplineBaseDeriv3;

}  // namespace base
}  // namespace sgpp

#endif /* LAGRANGE_SPLINE_BASE_DERIV3_HPP */
