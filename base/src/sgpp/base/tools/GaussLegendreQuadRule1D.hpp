// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GAUSSLEGENDREQUADRULE1D_HPP_
#define GAUSSLEGENDREQUADRULE1D_HPP_

#include <sgpp/base/tools/QuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

class GaussLegendreQuadRule1D : public QuadRule1D {
 public:
  /**
   * load gauss quadrature points for uniform weight function. The points
   * and the weights are generated with numpy.polynomial.legendre.leggauss.
   * the weights are additionally normalized to 1.
   */
  GaussLegendreQuadRule1D();
  ~GaussLegendreQuadRule1D() override;

  // delete the copy constructor
  GaussLegendreQuadRule1D(const GaussLegendreQuadRule1D& that) = delete;

  // get the maximum level that is supported by the quadrature rule
  size_t getMaxSupportedLevel() const override;

  /**
   * the coordinates are normalized to [0, 1].
   *
   * @param level level of quadrature, is equal to the number of quadrature points
   * @param coordinates returns the x-coordinates in [0, 1]
   * @param weights returns the corresponding weights (scaled by 0.5)
   */
  void getLevelPointsAndWeightsNormalized(size_t level, base::DataVector& coordinates,
                                          base::DataVector& weights);

  static GaussLegendreQuadRule1D& getInstance();

 protected:
  const size_t maxSupportedLevel = 101;

  inline void setGaussLegendreValuesForDegree0();
  inline void setGaussLegendreValuesForDegree1();
  inline void setGaussLegendreValuesForDegree2();
  inline void setGaussLegendreValuesForDegree3();
  inline void setGaussLegendreValuesForDegree4();
  inline void setGaussLegendreValuesForDegree5();
  inline void setGaussLegendreValuesForDegree6();
  inline void setGaussLegendreValuesForDegree7();
  inline void setGaussLegendreValuesForDegree8();
  inline void setGaussLegendreValuesForDegree9();
  inline void setGaussLegendreValuesForDegree10();
  inline void setGaussLegendreValuesForDegree11();
  inline void setGaussLegendreValuesForDegree12();
  inline void setGaussLegendreValuesForDegree13();
  inline void setGaussLegendreValuesForDegree14();
  inline void setGaussLegendreValuesForDegree15();
  inline void setGaussLegendreValuesForDegree16();
  inline void setGaussLegendreValuesForDegree17();
  inline void setGaussLegendreValuesForDegree18();
  inline void setGaussLegendreValuesForDegree19();
  inline void setGaussLegendreValuesForDegree20();
  inline void setGaussLegendreValuesForDegree21();
  inline void setGaussLegendreValuesForDegree22();
  inline void setGaussLegendreValuesForDegree23();
  inline void setGaussLegendreValuesForDegree24();
  inline void setGaussLegendreValuesForDegree25();
  inline void setGaussLegendreValuesForDegree26();
  inline void setGaussLegendreValuesForDegree27();
  inline void setGaussLegendreValuesForDegree28();
  inline void setGaussLegendreValuesForDegree29();
  inline void setGaussLegendreValuesForDegree30();
  inline void setGaussLegendreValuesForDegree31();
  inline void setGaussLegendreValuesForDegree32();
  inline void setGaussLegendreValuesForDegree33();
  inline void setGaussLegendreValuesForDegree34();
  inline void setGaussLegendreValuesForDegree35();
  inline void setGaussLegendreValuesForDegree36();
  inline void setGaussLegendreValuesForDegree37();
  inline void setGaussLegendreValuesForDegree38();
  inline void setGaussLegendreValuesForDegree39();
  inline void setGaussLegendreValuesForDegree40();
  inline void setGaussLegendreValuesForDegree41();
  inline void setGaussLegendreValuesForDegree42();
  inline void setGaussLegendreValuesForDegree43();
  inline void setGaussLegendreValuesForDegree44();
  inline void setGaussLegendreValuesForDegree45();
  inline void setGaussLegendreValuesForDegree46();
  inline void setGaussLegendreValuesForDegree47();
  inline void setGaussLegendreValuesForDegree48();
  inline void setGaussLegendreValuesForDegree49();
  inline void setGaussLegendreValuesForDegree50();
  inline void setGaussLegendreValuesForDegree51();
  inline void setGaussLegendreValuesForDegree52();
  inline void setGaussLegendreValuesForDegree53();
  inline void setGaussLegendreValuesForDegree54();
  inline void setGaussLegendreValuesForDegree55();
  inline void setGaussLegendreValuesForDegree56();
  inline void setGaussLegendreValuesForDegree57();
  inline void setGaussLegendreValuesForDegree58();
  inline void setGaussLegendreValuesForDegree59();
  inline void setGaussLegendreValuesForDegree60();
  inline void setGaussLegendreValuesForDegree61();
  inline void setGaussLegendreValuesForDegree62();
  inline void setGaussLegendreValuesForDegree63();
  inline void setGaussLegendreValuesForDegree64();
  inline void setGaussLegendreValuesForDegree65();
  inline void setGaussLegendreValuesForDegree66();
  inline void setGaussLegendreValuesForDegree67();
  inline void setGaussLegendreValuesForDegree68();
  inline void setGaussLegendreValuesForDegree69();
  inline void setGaussLegendreValuesForDegree70();
  inline void setGaussLegendreValuesForDegree71();
  inline void setGaussLegendreValuesForDegree72();
  inline void setGaussLegendreValuesForDegree73();
  inline void setGaussLegendreValuesForDegree74();
  inline void setGaussLegendreValuesForDegree75();
  inline void setGaussLegendreValuesForDegree76();
  inline void setGaussLegendreValuesForDegree77();
  inline void setGaussLegendreValuesForDegree78();
  inline void setGaussLegendreValuesForDegree79();
  inline void setGaussLegendreValuesForDegree80();
  inline void setGaussLegendreValuesForDegree81();
  inline void setGaussLegendreValuesForDegree82();
  inline void setGaussLegendreValuesForDegree83();
  inline void setGaussLegendreValuesForDegree84();
  inline void setGaussLegendreValuesForDegree85();
  inline void setGaussLegendreValuesForDegree86();
  inline void setGaussLegendreValuesForDegree87();
  inline void setGaussLegendreValuesForDegree88();
  inline void setGaussLegendreValuesForDegree89();
  inline void setGaussLegendreValuesForDegree90();
  inline void setGaussLegendreValuesForDegree91();
  inline void setGaussLegendreValuesForDegree92();
  inline void setGaussLegendreValuesForDegree93();
  inline void setGaussLegendreValuesForDegree94();
  inline void setGaussLegendreValuesForDegree95();
  inline void setGaussLegendreValuesForDegree96();
  inline void setGaussLegendreValuesForDegree97();
  inline void setGaussLegendreValuesForDegree98();
  inline void setGaussLegendreValuesForDegree99();
  inline void setGaussLegendreValuesForDegree100();
};

}  // namespace base
}  // namespace sgpp

#endif /* GAUSSLEGENDREQUADRULE1D_HPP_ */
