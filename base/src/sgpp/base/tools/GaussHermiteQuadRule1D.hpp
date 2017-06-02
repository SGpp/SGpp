// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GAUSSHERMITEQUADRULE1D_HPP_
#define GAUSSHERMITEQUADRULE1D_HPP_

#include <sgpp/base/tools/QuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

class GaussHermiteQuadRule1D : public QuadRule1D {
 public:
  GaussHermiteQuadRule1D();
  ~GaussHermiteQuadRule1D() override;

  // delete the copy constructor
  GaussHermiteQuadRule1D(const GaussHermiteQuadRule1D& that) = delete;

  /**
   * load gauss quadrature points for standard normal weight function. The points
   * and the weights are generated with numpy.polynomial.hermite.hermgauss,
   * the coordinates are scaled by sqrt(2), the weights are normalized to 1.
   */

  /**
   * the coordinates are scaled by sqrt(2) and then normalized with respect
   * to a given mean and standard deviation. The weights are normalized
   * to 1.
   *
   * @param level level of quadrature, is equal to the number of quadrature points
   * @param coordinates returns the x-coordinates in [-infty, infty]
   * @param weights returns the corresponding weights (scaled by sqrt(2))
   * @param mean mean of the normal distribution the coordinates should be transformed to
   * @param stdd standard deviation of the normal distribution the coordinates should be transformed
   * to
   */
  void getLevelPointsAndWeightsNormalized(size_t level, base::DataVector& coordinates,
                                          base::DataVector& weights, double mean = 0.0f,
                                          double stdd = 1.0f);

  static GaussHermiteQuadRule1D& getInstance();

 protected:
  inline void setGaussHermiteValuesForDegree0();
  inline void setGaussHermiteValuesForDegree1();
  inline void setGaussHermiteValuesForDegree2();
  inline void setGaussHermiteValuesForDegree3();
  inline void setGaussHermiteValuesForDegree4();
  inline void setGaussHermiteValuesForDegree5();
  inline void setGaussHermiteValuesForDegree6();
  inline void setGaussHermiteValuesForDegree7();
  inline void setGaussHermiteValuesForDegree8();
  inline void setGaussHermiteValuesForDegree9();
  inline void setGaussHermiteValuesForDegree10();
  inline void setGaussHermiteValuesForDegree11();
  inline void setGaussHermiteValuesForDegree12();
  inline void setGaussHermiteValuesForDegree13();
  inline void setGaussHermiteValuesForDegree14();
  inline void setGaussHermiteValuesForDegree15();
  inline void setGaussHermiteValuesForDegree16();
  inline void setGaussHermiteValuesForDegree17();
  inline void setGaussHermiteValuesForDegree18();
  inline void setGaussHermiteValuesForDegree19();
  inline void setGaussHermiteValuesForDegree20();
  inline void setGaussHermiteValuesForDegree21();
  inline void setGaussHermiteValuesForDegree22();
  inline void setGaussHermiteValuesForDegree23();
  inline void setGaussHermiteValuesForDegree24();
  inline void setGaussHermiteValuesForDegree25();
  inline void setGaussHermiteValuesForDegree26();
  inline void setGaussHermiteValuesForDegree27();
  inline void setGaussHermiteValuesForDegree28();
  inline void setGaussHermiteValuesForDegree29();
  inline void setGaussHermiteValuesForDegree30();
  inline void setGaussHermiteValuesForDegree31();
  inline void setGaussHermiteValuesForDegree32();
  inline void setGaussHermiteValuesForDegree33();
  inline void setGaussHermiteValuesForDegree34();
  inline void setGaussHermiteValuesForDegree35();
  inline void setGaussHermiteValuesForDegree36();
  inline void setGaussHermiteValuesForDegree37();
  inline void setGaussHermiteValuesForDegree38();
  inline void setGaussHermiteValuesForDegree39();
  inline void setGaussHermiteValuesForDegree40();
  inline void setGaussHermiteValuesForDegree41();
  inline void setGaussHermiteValuesForDegree42();
  inline void setGaussHermiteValuesForDegree43();
  inline void setGaussHermiteValuesForDegree44();
  inline void setGaussHermiteValuesForDegree45();
  inline void setGaussHermiteValuesForDegree46();
  inline void setGaussHermiteValuesForDegree47();
  inline void setGaussHermiteValuesForDegree48();
  inline void setGaussHermiteValuesForDegree49();
};

}  // namespace base
}  // namespace sgpp

#endif /* GAUSSHERMITEQUADRULE1D_HPP_ */
