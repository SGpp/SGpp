/*
 * AnisotropyLoadModel.hpp
 *
 *  Created on: Oct 9, 2013
 *      Author: heenemo
 */

#ifndef ANISOTROPYLOADMODEL_HPP_
#define ANISOTROPYLOADMODEL_HPP_

#include <stddef.h>
#include <cmath>
#include <vector>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"

namespace combigrid {

class AnisotropyLoadModel: public LoadModel {
 public:
  AnisotropyLoadModel();

  inline real eval(const LevelVector& l) const;

  virtual ~AnisotropyLoadModel();
};

inline real AnisotropyLoadModel::eval(const LevelVector& l) const {
  // number of grid points
  LevelType lsum = sum(l);

  // todo: remove 0.5 - war nur wegen konsistenz mit matlab implementierung
  const real N = std::pow(real(2.0), static_cast<real>(lsum));

  // coefficients for iso model 1.252e-05,1.076,0.01308
  const real m(1.252e-05);
  const real k(1.076);
  const real ci(0.01308);

  // value of iso model
  const real r = m * std::pow(N, k) + ci;

  // coefficients for p=2 aniso model
  const real c11(15.3450);
  const real c21(1.0204);
  const real c31(2.9396);
  const real c1(-10.7883);
  const real c22(2.3668);
  const real c32(-1.1119);
  const real c2(-0.8959);
  const real c33(0.3434);
  const real c3(-0.4814);
  const real c(2.6451);

  // anisotropy vector
  std::vector<real> s(l.size(), 0.0);

  for (size_t i = 0; i < s.size(); ++i)
    s[i] = static_cast<real>(l[i]) / static_cast<real>(lsum);

  // variables for aniso model
  const real s1 = s[0];
  const real s2 = s[2];
  const real s3 = s[3];

  // constant contribution
  real h = c;

  // linear contribution
  h += c1 * s1 + c2 * s2 + c3 * s3;

  // p=2 polynomial contribution
  h += c11 * s1 * s1 + c21 * s2 * s1 + c31 * s3 * s1 + c22 * s2 * s2
       + c32 * s3 * s2 + c33 * s3 * s3;

  return r * h;
}

} /* namespace combigrid */
#endif /* ANISOTROPYLOADMODEL_HPP_ */
