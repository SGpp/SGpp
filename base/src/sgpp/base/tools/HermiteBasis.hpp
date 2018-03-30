// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HERMITEBASIS_HPP
#define HERMITEBASIS_HPP

namespace sgpp {
namespace base {
class HermiteBasis {
  // hermite basis functions according to
  // https://en.wikipedia.org/wiki/Cubic_Hermite_spline
 public:
  static double h_0_0(double x) { return (1+2*x)*(1-x)*(1-x);}
  static double h_1_0(double x) { return x*(1-x)*(1-x);}
  static double h_0_1(double x) { return x*x*(3 - 2*x);}
  static double h_1_1(double x) { return x*x*(x-1);}
};

}  // namespace base
}  // namespace sgpp

#endif /* HERMITEBASIS_HPP */
