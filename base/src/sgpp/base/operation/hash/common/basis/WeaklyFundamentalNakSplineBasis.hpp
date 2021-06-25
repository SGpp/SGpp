// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_HPP
#define WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental not-a-knot spline basis.
 */
template <class LT, class IT>
class WeaklyFundamentalNakSplineBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalNakSplineBasis()
      : weaklyFundamentalSplineBasis(WeaklyFundamentalSplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalNakSplineBasis(size_t degree)
      : weaklyFundamentalSplineBasis(WeaklyFundamentalSplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~WeaklyFundamentalNakSplineBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    double t = x * static_cast<double>(hInv) - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        return weaklyFundamentalSplineBasis.eval(l, i, x);

      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return 1.0 - t * t;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return weaklyFundamentalSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 2.0192307692307693e-01;
              result = -9.5192307692307687e-01 + result * t;
              result = 1.0961538461538463e+00 + result * t;
              result *= t;
              return result;
            } else {
              t -= 1.0;
              double result = -5.7692307692307696e-02;
              result = 2.5961538461538464e-01 + result * t;
              result = -2.8846153846153844e-01 + result * t;
              result *= t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 2.1428571428571427e-01;
              result = -9.6428571428571430e-01 + result * t;
              result = 1.0714285714285714e+00 + result * t;
              result *= t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.2500000000000000e-01;
              result = 3.2142857142857145e-01 + result * t;
              result = -2.1428571428571427e-01 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = 1.7857142857142856e-02;
              result = -5.3571428571428568e-02 + result * t;
              result = 5.3571428571428568e-02 + result * t;
              result = -1.7857142857142856e-02 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 1.2500000000000000e-01;
              result = -1.8750000000000000e-01 + result * t;
              result = -1.2500000000000000e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -5.9375000000000000e-01;
              result = 5.6250000000000000e-01 + result * t;
              result = 6.2500000000000000e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 6.5625000000000000e-01;
              result = -1.2187500000000000e+00 + result * t;
              result = -3.1250000000000000e-02 + result * t;
              result = 5.9375000000000000e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -2.9166666666666669e-01;
              result = 7.5000000000000000e-01 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = 4.1666666666666664e-02;
              result = -1.2500000000000000e-01 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              result = -4.1666666666666664e-02 + result * t;
              return result;
            }
          }
        }

      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return 1.0 - t * t;
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return weaklyFundamentalSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -1.6666666666666666e-01;
            result = 8.3333333333333337e-01 + result * t;
            result = -8.3333333333333337e-01 + result * t;
            result = -8.3333333333333337e-01 + result * t;
            result = 1.0000000000000000e+00 + result * t;
            return result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 4.1145470470972381e-03;
              result = -5.2880937409913717e-02 + result * t;
              result = 2.5243387314625304e-01 + result * t;
              result = -5.2691055961225519e-01 + result * t;
              result = 4.0130037316525197e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -2.0620724676844440e-03;
              result = 8.8372682965448569e-03 + result * t;
              result = -1.1828141533960099e-02 + result * t;
              result = 3.5138128493598849e-04 + result * t;
              result = 1.0801904244253702e-02 + result * t;
              result = -6.1003398240900043e-03 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 2.9583891062696650e-04;
              result = -1.4730940418773624e-03 + result * t;
              result = 2.9002069753748905e-03 + result * t;
              result = -2.7301582145196038e-03 + result * t;
              result = 1.0589530600025915e-03 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -4.3062432405671986e-07;
              result = 6.1005112574701972e-06 + result * t;
              result = -3.3780085864893800e-05 + result * t;
              result = 9.0287566610558926e-05 + result * t;
              result = -1.1392405728656110e-04 + result * t;
              result = 5.1746689607482503e-05 + result * t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.3884771234422369e-02;
              result = -1.1970944282074093e-01 + result * t;
              result = 2.5338631824956254e-01 + result * t;
              result = 1.7139797018949438e-01 + result * t;
              result = -6.2082201056206943e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = -3.4307459460041163e-02;
              result = 8.8562125695594601e-02 + result * t;
              result = 6.6502415498684575e-02 + result * t;
              result = -2.6354684459041328e-01 + result * t;
              result = -5.6291071385875184e-02 + result * t;
              result = 1.9908083424205042e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.9772705970710753e-02;
              result = -8.2975171604611195e-02 + result * t;
              result = 7.7676323680651374e-02 + result * t;
              result = 1.2425856147879648e-01 + result * t;
              result = -2.0116630858847537e-01 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.4633576017002452e-03;
              result = 1.5888358248942556e-02 + result * t;
              result = -5.6497303030685918e-02 + result * t;
              result = 5.7163562600190905e-02 + result * t;
              result = 4.7342628846180658e-02 + result * t;
              result = -6.2433889062927968e-02 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 4.1250937315480164e-03;
              result = -5.2940966935253082e-02 + result * t;
              result = 2.5231949105027229e-01 + result * t;
              result = -5.2580228771389204e-01 + result * t;
              result = 3.9985284700395129e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -2.1055637597912490e-03;
              result = 8.9354390379671640e-03 + result * t;
              result = -1.1713676333443215e-02 + result * t;
              result = 3.6224752856623638e-05 + result * t;
              result = 1.0703911347564831e-02 + result * t;
              result = -5.8563350451541546e-03 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 3.2909181727107666e-04;
              result = -1.5923797609890807e-03 + result * t;
              result = 2.9724422205129504e-03 + result * t;
              result = -2.5478076175825290e-03 + result * t;
              result = 8.4926920586084306e-04 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -1.0615865073260537e-05;
              result = 5.3079325366302692e-05 + result * t;
              result = -1.0615865073260538e-04 + result * t;
              result = 1.0615865073260538e-04 + result * t;
              result = -5.3079325366302692e-05 + result * t;
              result = 1.0615865073260537e-05 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.3799337496874365e-02;
              result = -1.1764498958355420e-01 + result * t;
              result = 2.4592473021998218e-01 + result * t;
              result = 1.6902170570373345e-01 + result * t;
              result = -6.0137181556895192e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = -3.7647688975389365e-02;
              result = 8.9345072869561282e-02 + result * t;
              result = 7.6125229936024733e-02 + result * t;
              result = -2.4466403567227474e-01 + result * t;
              result = -6.4201054196767235e-02 + result * t;
              result = 1.8104247603884530e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 2.8217942728735564e-02;
              result = -9.8893372007385516e-02 + result * t;
              result = 5.7028631660376264e-02 + result * t;
              result = 1.4330520159927357e-01 + result * t;
              result = -1.5601158913194416e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -9.8889538888502406e-03;
              result = 4.2196341636292302e-02 + result * t;
              result = -5.6365429081810178e-02 + result * t;
              result = 3.2102918234448797e-03 + result * t;
              result = 4.7200934661867512e-02 + result * t;
              result = -2.6353185150944279e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.4980084136448403e-03;
              result = -7.2484278079589038e-03 + result * t;
              result = 1.3530398574856621e-02 + result * t;
              result = -1.1597484492734246e-02 + result * t;
              result = 3.8658281642447489e-03 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -4.8322852053059363e-05;
              result = 2.4161426026529680e-04 + result * t;
              result = -4.8322852053059361e-04 + result * t;
              result = 4.8322852053059361e-04 + result * t;
              result = -2.4161426026529680e-04 + result * t;
              result = 4.8322852053059363e-05 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 3.6044172613134747e-03;
              result = -1.4909145287133148e-02 + result * t;
              result = -6.0604127900827486e-03 + result * t;
              result = 2.8393347730744867e-02 + result * t;
              result = 2.9057441814890835e-02 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -5.3234577203784571e-02;
              result = 3.9157113632568977e-02 + result * t;
              result = 1.3942739728253223e-01 + result * t;
              result = 1.4194844766944842e-01 + result * t;
              result = -1.1461231731129669e-01 + result * t;
              result = -1.5268606406946839e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 1.1753404783757340e-01;
              result = -2.2701577238635387e-01 + result * t;
              result = -2.3628992022503756e-01 + result * t;
              result = 2.6282754927461327e-01 + result * t;
              result = 4.7802233838654995e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = -1.2833863685436561e-01;
              result = 3.6065446680151309e-01 + result * t;
              result = 3.0987468605280882e-02 + result * t;
              result = -6.3279636734288869e-01 + result * t;
              result = -2.5585174096884770e-02 + result * t;
              result = 3.9507824288734517e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 7.8473853438934288e-02;
              result = -2.8103871747031500e-01 + result * t;
              result = 1.9021896726767709e-01 + result * t;
              result = 3.4070647073837634e-01 + result * t;
              result = -3.9729082003259530e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -2.6065958111654366e-02;
              result = 1.1133054972435644e-01 + result * t;
              result = -1.4919736822424004e-01 + result * t;
              result = 9.8696021088604830e-03 + result * t;
              result = 1.2299342056060007e-01 + result * t;
              result = -6.8930246057922598e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 3.9265097723425145e-03;
              result = -1.8999240833915392e-02 + result * t;
              result = 3.5465249556642063e-02 + result * t;
              result = -3.0398785334264627e-02 + result * t;
              result = 1.0132928444754876e-02 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -1.2666160555943596e-04;
              result = 6.3330802779717978e-04 + result * t;
              result = -1.2666160555943596e-03 + result * t;
              result = 1.2666160555943596e-03 + result * t;
              result = -6.3330802779717978e-04 + result * t;
              result = 1.2666160555943596e-04 + result * t;
              return result;
            }
          }
        }

      case 7:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return 1.0 - t * t;
        } else if ((i > 9) && (i < hInv - 9)) {
          // l >= 5, 9 < i < 2^l - 9
          return weaklyFundamentalSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -1.6666666666666666e-01;
            result = 8.3333333333333337e-01 + result * t;
            result = -8.3333333333333337e-01 + result * t;
            result = -8.3333333333333337e-01 + result * t;
            result = 1.0000000000000000e+00 + result * t;
            return result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 2.2503649501781186e-05;
              result = -6.0667311908106972e-04 + result * t;
              result = 6.8269755058389522e-03 + result * t;
              result = -4.0835982396216351e-02 + result * t;
              result = 1.3582445212555574e-01 + result * t;
              result = -2.3523375669133104e-01 + result * t;
              result = 1.6259926219922699e-01 + result * t;
              result *= t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.1754318967321538e-06;
              result = 2.3429066968803494e-05 + result * t;
              result = -1.7195311950824255e-04 + result * t;
              result = 5.1015402509581493e-04 + result * t;
              result = -1.4407816748338684e-04 + result * t;
              result = -2.0366158695079022e-03 + result * t;
              result = 2.6450596187038589e-03 + result * t;
              result *= t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 3.9587351846363302e-05;
              result = -8.7259754269615290e-04 + result * t;
              result = 7.0067813511653382e-03 + result * t;
              result = -2.1711347330167205e-02 + result * t;
              result = -2.3700033078373774e-03 + result * t;
              result = 1.3432426280688878e-01 + result * t;
              result = -1.7219670451162619e-01 + result * t;
              result *= t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.8319997233243217e-05;
              result = 2.3584830900201965e-04 + result * t;
              result = -6.3420945316426029e-04 + result * t;
              result = -2.3234624180833233e-03 + result * t;
              result = 9.1112734882810315e-03 + result * t;
              result = 6.4467943135746832e-03 + result * t;
              result = -2.6978298050115802e-02 + result * t;
              result *= t;
              return result;
            }
          } else if ((l == 4) && (i == 7)) {
            // l = 4, i = 7
            if ((t < -7.0) || (t > 9.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 4.4756604140043971e-05;
              result = -3.5474922508209859e-04 + result * t;
              result = 1.6336284038032301e-04 + result * t;
              result = 2.2623931515483426e-03 + result * t;
              result = 1.3813268802489874e-03 + result * t;
              result = -5.5122086760940213e-03 + result * t;
              result = -6.7262882896154744e-03 + result * t;
              result *= t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -1.8932333752693443e-03;
              result = 8.9843569083913259e-04 + result * t;
              result = 6.6876004294647304e-03 + result * t;
              result = 2.0644629213149631e-02 + result * t;
              result = 1.0657836755581922e-02 + result * t;
              result = -6.6985334608811620e-02 + result * t;
              result = -9.2561069964300122e-02 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 6.4242576114280758e-03;
              result = -1.2354197936046278e-02 + result * t;
              result = -2.7679686306156703e-02 + result * t;
              result = 1.2959985886332251e-03 + result * t;
              result = 1.1181790358518336e-01 + result * t;
              result = 1.2945058971341000e-01 + result * t;
              result = -8.6403729397106027e-02 + result * t;
              result = -1.2255113585934567e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.1137563913351757e-02;
              result = 3.2615605343950256e-02 + result * t;
              result = 3.3104535917555231e-02 + result * t;
              result = -9.7566385582861778e-02 + result * t;
              result = -1.8202990744279365e-01 + result * t;
              result = 1.4547986973848787e-01 + result * t;
              result = 3.4558133927273232e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 1.1945996989810794e-02;
              result = -4.5347342049512039e-02 + result * t;
              result = -5.0906741991301073e-03 + result * t;
              result = 1.6737563719685677e-01 + result * t;
              result = 2.1247279313005213e-02 + result * t;
              result = -3.9961756893264450e-01 + result * t;
              result = -1.6560821652104624e-02 + result * t;
              result = 2.6604749333371847e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -8.3163814615979416e-03;
              result = 3.8274636879163527e-02 + result * t;
              result = -2.6308789710175629e-02 + result * t;
              result = -1.2017796989809650e-01 + result * t;
              result = 1.5100613976226829e-01 + result * t;
              result = 1.8812715623955681e-01 + result * t;
              result = -2.9646701715499812e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 3.7087657052593680e-03;
              result = -1.9940033352022060e-02 + result * t;
              result = 2.8695020871248776e-02 + result * t;
              result = 3.1324283582550350e-02 + result * t;
              result = -1.1837425050453133e-01 + result * t;
              result = 5.6465401529922604e-02 + result * t;
              result = 9.1983037511451843e-02 + result * t;
              result = -7.3862225343879540e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -9.8041556667037891e-04;
              result = 6.0213265847935144e-03 + result * t;
              result = -1.3061099430436859e-02 + result * t;
              result = 5.5056873425412042e-03 + result * t;
              result = 2.4879225181794497e-02 + result * t;
              result = -4.4977860245765722e-02 + result * t;
              result = 2.4884487568831534e-02 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.2072422471724799e-04;
              result = -8.4158238189913710e-04 + result * t;
              result = 2.4781331782462732e-03 + result * t;
              result = -3.7944558712036326e-03 + result * t;
              result = 2.4029671099977529e-03 + result * t;
              result = 1.8141169223211654e-03 + result * t;
              result = -4.4512546172674602e-03 + result * t;
              result = 2.2713514350877913e-03 + result * t;
              return result;
            } else {
              t -= 5.0;
              double result = -1.3914756881238442e-07;
              result = 3.4871911215987656e-06 + result * t;
              result = -3.6152394086341924e-05 + result * t;
              result = 1.9782215664435564e-04 + result * t;
              result = -5.9982436523310959e-04 + result * t;
              result = 9.4908779813050975e-04 + result * t;
              result = -6.0570175458912344e-04 + result * t;
              result *= t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 2.6036925510092182e-05;
              result = -6.5237050374241724e-04 + result * t;
              result = 6.7615644427544194e-03 + result * t;
              result = -3.6988327885346672e-02 + result * t;
              result = 1.1212038017839150e-01 + result * t;
              result = -1.7735150538327948e-01 + result * t;
              result = 1.1315257493880702e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -1.3067934122423832e-05;
              result = 7.6663410540163884e-05 + result * t;
              result = -1.4692067567262032e-04 + result * t;
              result = -3.2467858319237330e-06 + result * t;
              result = 4.1405263368380322e-04 + result * t;
              result = -5.8986507489320434e-04 + result * t;
              result = 2.7568174422425856e-04 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 2.2552954658425583e-06;
              result = -1.4812128316802933e-05 + result * t;
              result = 3.8633170997462535e-05 + result * t;
              result = -4.5276700377401174e-05 + result * t;
              result = 7.7492501483486359e-06 + result * t;
              result = 3.9129895972017469e-05 + result * t;
              result = -4.0976101817520529e-05 + result * t;
              result = 1.3297317928053443e-05 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -1.4038250895698593e-07;
              result = 9.7493994409497317e-07 + result * t;
              result = -2.8783941206613495e-06 + result * t;
              result = 4.6425711623570155e-06 + result * t;
              result = -4.3330664181998812e-06 + result * t;
              result = 2.2284341579313672e-06 + result * t;
              result = -4.9520759065141497e-07 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 1.1053740862754798e-09;
              result = -7.7376186039283589e-09 + result * t;
              result = 2.3212855811785075e-08 + result * t;
              result = -3.8688093019641789e-08 + result * t;
              result = 3.8688093019641789e-08 + result * t;
              result = -2.3212855811785075e-08 + result * t;
              result = 7.7376186039283589e-09 + result * t;
              result = -1.1053740862754798e-09 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 2.2108026631303853e-04;
              result = -4.3130889168398112e-03 + result * t;
              result = 3.0385702223110777e-02 + result * t;
              result = -7.9991415288870771e-02 + result * t;
              result = -3.6061653696913241e-02 + result * t;
              result = 5.0343956936202849e-01 + result * t;
              result = -5.8500272890037075e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -5.0163220710328841e-04;
              result = 1.8771585399252673e-03 + result * t;
              result = 1.1545377001362510e-03 + result * t;
              result = -1.2198914327003638e-02 + result * t;
              result = 5.9134299887455260e-03 + result * t;
              result = 3.0221886189078145e-02 + result * t;
              result = -3.3353929161460835e-02 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 3.1944023124356557e-04;
              result = -1.6342669097977518e-03 + result * t;
              result = 1.8832125905187979e-03 + result * t;
              result = 4.1740250239415324e-03 + result * t;
              result = -1.1350806768016261e-02 + result * t;
              result = 3.9371689443653641e-03 + result * t;
              result = 9.5586901654273245e-03 + result * t;
              result = -6.8874632776825713e-03 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -1.0103098001508067e-04;
              result = 6.0181470890720735e-04 + result * t;
              result = -1.2141440121528357e-03 + result * t;
              result = 2.5649242309404058e-04 + result * t;
              result = 2.6724891305076079e-03 + result * t;
              result = -4.0447341016976446e-03 + result * t;
              result = 1.9232509583878375e-03 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.6044157194106501e-05;
              result = -1.0540215119835737e-04 + result * t;
              result = 2.7509366097371415e-04 + result * t;
              result = -3.2309130458985140e-04 + result * t;
              result = 5.7228578971736657e-05 + result * t;
              result = 2.7581776015248239e-04 + result * t;
              result = -2.8982882853496373e-04 + result * t;
              result = 9.4138127031132804e-05 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -9.9453803103094462e-07;
              result = 6.9069491603881354e-06 + result * t;
              result = -2.0391945140193545e-05 + result * t;
              result = 3.2890234097086359e-05 + result * t;
              result = -3.0697551823947269e-05 + result * t;
              result = 1.5787312366601453e-05 + result * t;
              result = -3.5082916370225450e-06 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 7.8310081183538954e-09;
              result = -5.4817056828477266e-08 + result * t;
              result = 1.6445117048543181e-07 + result * t;
              result = -2.7408528414238632e-07 + result * t;
              result = 2.7408528414238632e-07 + result * t;
              result = -1.6445117048543181e-07 + result * t;
              result = 5.4817056828477266e-08 + result * t;
              result = -7.8310081183538954e-09 + result * t;
              return result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 2.6216227411080255e-04;
              result = -3.5998840686331641e-03 + result * t;
              result = 1.3199581657453090e-02 + result * t;
              result = 9.0533481228980296e-03 + result * t;
              result = -9.2870661409871463e-02 + result * t;
              result = -2.3196224441615544e-02 + result * t;
              result = 2.3267290767345314e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -2.1756203320584594e-03;
              result = 3.7406596064693078e-03 + result * t;
              result = 1.4888888111486814e-02 + result * t;
              result = -3.6837011918017877e-03 + result * t;
              result = -9.4961658068667479e-02 + result * t;
              result = -6.8077618645357770e-03 + result * t;
              result = 2.0146201109132153e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 3.3966566834791193e-03;
              result = -1.1488682717939909e-02 + result * t;
              result = -8.3551812229249872e-03 + result * t;
              result = 5.0723921840625817e-02 + result * t;
              result = 3.7858898786333584e-02 + result * t;
              result = -1.5448419498266885e-01 + result * t;
              result = -3.0114235739118922e-02 + result * t;
              result = 1.1246281735221415e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -2.7724353053418401e-03;
              result = 1.2287914066413925e-02 + result * t;
              result = -5.9574871775029351e-03 + result * t;
              result = -4.4499241121328574e-02 + result * t;
              result = 4.6312103482557987e-02 + result * t;
              result = 7.8883769774799814e-02 + result * t;
              result = -1.0954164762086316e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.3516194313367197e-03;
              result = -7.1191330709789551e-03 + result * t;
              result = 9.5488558088019770e-03 + result * t;
              result = 1.2996798300401232e-02 + result * t;
              result = -4.2536687136471546e-02 + result * t;
              result = 1.7347331303503220e-02 + result * t;
              result = 3.3698239264672129e-02 + result * t;
              result = -2.5287023901264777e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -3.9175625128765431e-04;
              result = 2.3422029483780838e-03 + result * t;
              result = -4.7819345590006355e-03 + result * t;
              result = 1.2607613765119885e-03 + result * t;
              result = 9.8630828303592454e-03 + result * t;
              result = -1.5196370222097471e-02 + result * t;
              result = 7.2608503013620491e-03 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 6.0896723955089969e-05;
              result = -4.0009081063549658e-04 + result * t;
              result = 1.0444018542271256e-03 + result * t;
              result = -1.2273359878878364e-03 + result * t;
              result = 2.1937291889461460e-04 + result * t;
              result = 1.0442638866763549e-03 + result * t;
              result = -1.0983450094554580e-03 + result * t;
              result = 3.5683642422560611e-04 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -3.7705834981484294e-06;
              result = 2.6186257050133186e-05 + result * t;
              result = -7.7311806528964646e-05 + result * t;
              result = 1.2469646214349138e-04 + result * t;
              result = -1.1638336466725860e-04 + result * t;
              result = 5.9854301828875857e-05 + result * t;
              result = -1.3300955961972412e-05 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 2.9689633843688420e-08;
              result = -2.0782743690581893e-07 + result * t;
              result = 6.2348231071745685e-07 + result * t;
              result = -1.0391371845290946e-06 + result * t;
              result = 1.0391371845290946e-06 + result * t;
              result = -6.2348231071745685e-07 + result * t;
              result = 2.0782743690581893e-07 + result * t;
              result = -2.9689633843688420e-08 + result * t;
              return result;
            }
          } else if (i == 7) {
            // l >= 5, i = 7
            if ((t < -7.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 4.4735023389159415e-05;
              result = -3.5456764543243316e-04 + result * t;
              result = 1.6326443418213484e-04 + result * t;
              result = 2.2611957232376006e-03 + result * t;
              result = 1.3806292197650139e-03 + result * t;
              result = -5.5092426891575974e-03 + result * t;
              result = -6.7227050766281628e-03 + result * t;
              result *= t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -1.8940626123998068e-03;
              result = 8.9801300946403039e-04 + result * t;
              result = 6.6846088025613011e-03 + result * t;
              result = 2.0636701894813418e-02 + result * t;
              result = 1.0661293674062058e-02 + result * t;
              result = -6.6935480244660894e-02 + result * t;
              result = -9.2502886062100734e-02 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 6.4357508804497726e-03;
              result = -1.2360425277334617e-02 + result * t;
              result = -2.7702628001050458e-02 + result * t;
              result = 1.2377496155871410e-03 + result * t;
              result = 1.1172225803421611e-01 + result * t;
              result = 1.2940958045358331e-01 + result * t;
              result = -8.6290474167190628e-02 + result * t;
              result = -1.2245181153826062e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.1188285488039624e-02;
              result = 3.2689830885813786e-02 + result * t;
              result = 3.3285588824387057e-02 + result * t;
              result = -9.7430488733942375e-02 + result * t;
              result = -1.8231024824489023e-01 + result * t;
              result = 1.4472096156867587e-01 + result * t;
              result = 3.4502102379886129e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 1.2068149738337627e-02;
              result = -4.5628167530463573e-02 + result * t;
              result = -5.5294211095622998e-03 + result * t;
              result = 1.6775492659381294e-01 + result * t;
              result = 2.3030310698099812e-02 + result * t;
              result = -3.9854335928740375e-01 + result * t;
              result = -1.7940821713686512e-02 + result * t;
              result = 2.6478838261086579e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -8.5013647890198798e-03;
              result = 3.8848880637899812e-02 + result * t;
              result = -2.5867281787253591e-02 + result * t;
              result = -1.2192945106913526e-01 + result * t;
              result = 1.4857769621027397e-01 + result * t;
              result = 1.9079155282228680e-01 + result * t;
              result = -2.9185596438117245e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 3.8933494692219724e-03;
              result = -2.0660672885239349e-02 + result * t;
              result = 2.8697341470727791e-02 + result * t;
              result = 3.3919581947398143e-02 + result * t;
              result = -1.1838308079650259e-01 + result * t;
              result = 5.0479666164840870e-02 + result * t;
              result = 9.1989746985673748e-02 + result * t;
              result = -6.9935932356120584e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -1.1015228993508576e-03;
              result = 6.5927733993144553e-03 + result * t;
              result = -1.3506356987046893e-02 + result * t;
              result = 3.7634274452158743e-03 + result * t;
              result = 2.7322435418349913e-02 + result * t;
              result = -4.2328424257928980e-02 + result * t;
              result = 2.0254281042196943e-02 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.7014672400308872e-04;
              result = -1.1178868961415481e-03 + result * t;
              result = 2.9183025224718290e-03 + result * t;
              result = -3.4300579775817766e-03 + result * t;
              result = 6.1474183775357140e-04 + result * t;
              result = 2.9154969012958956e-03 + result * t;
              result = -3.0673562725515159e-03 + result * t;
              result = 9.9661316075045627e-04 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -1.0531521347810959e-05;
              result = 7.3140171880072967e-05 + result * t;
              result = -2.1593765031259636e-04 + result * t;
              result = 3.4828653276225218e-04 + result * t;
              result = -3.2506743057810205e-04 + result * t;
              result = 1.6717753572588107e-04 + result * t;
              result = -3.7150563494640233e-05 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 8.2925364943393381e-08;
              result = -5.8047755460375364e-07 + result * t;
              result = 1.7414326638112610e-06 + result * t;
              result = -2.9023877730187685e-06 + result * t;
              result = 2.9023877730187685e-06 + result * t;
              result = -1.7414326638112610e-06 + result * t;
              result = 5.8047755460375364e-07 + result * t;
              result = -8.2925364943393381e-08 + result * t;
              return result;
            }
          } else {
            // l >= 5, i = 9
            if ((t < -9.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -5.0) {
              t += 9.0;
              double result = 1.3780287234530888e-07;
              result = -4.0498887655435551e-07 + result * t;
              result = -7.7900805975229674e-07 + result * t;
              result = 9.4972858589513266e-09 + result * t;
              result = 2.1891318473151508e-06 + result * t;
              result = 3.1441442051516831e-06 + result * t;
              result = 1.4835950892407763e-06 + result * t;
              result *= t;
              return result;
            } else if (t < -4.0) {
              t += 5.0;
              double result = -1.1957354658434113e-04;
              result = 3.4534915491142933e-06 + result * t;
              result = 3.5803024010966957e-05 + result * t;
              result = 1.9591043977125960e-04 + result * t;
              result = 5.9402777308508348e-04 + result * t;
              result = 9.3991598851872005e-04 + result * t;
              result = 5.9984836443326931e-04 + result * t;
              result *= t;
              return result;
            } else if (t < -3.0) {
              t += 4.0;
              double result = 9.7297028475482625e-04;
              result = -8.3356133454127361e-04 + result * t;
              result = -2.4545205049655112e-03 + result * t;
              result = -3.7583461973891307e-03 + result * t;
              result = -2.3803045271898623e-03 + result * t;
              result = 1.7962500814767482e-03 + result * t;
              result = 4.4081266630701310e-03 + result * t;
              result = 2.2493855347840725e-03 + result * t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -3.7057082860908006e-03;
              result = 5.9772306587425102e-03 + result * t;
              result = 1.2976487467638197e-02 + result * t;
              result = 5.5195912260831270e-03 + result * t;
              result = -2.4576161090808051e-02 + result * t;
              result = -4.4510989772350490e-02 + result * t;
              result = -2.4636850083893897e-02 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 8.4480795759128233e-03;
              result = -1.9962727343893094e-02 + result * t;
              result = -2.8980002587813553e-02 + result * t;
              result = 3.0360698432233747e-02 + result * t;
              result = 1.1711190165157861e-01 + result * t;
              result = 5.6481534861336939e-02 + result * t;
              result = -9.0503084708676074e-02 + result * t;
              result = -7.2956399880679410e-02 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.2577812087959618e-02;
              result = 3.9173829687496664e-02 + result * t;
              result = 2.8653304442997164e-02 + result * t;
              result = -1.1829743950828163e-01 + result * t;
              result = -1.5481709221853499e-01 + result * t;
              result = 1.7815016546711263e-01 + result * t;
              result = 2.8969866372663206e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 1.2649073554246904e-02;
              result = -4.8870854928220658e-02 + result * t;
              result = -4.3777127917480577e-04 + result * t;
              result = 1.7235310494056758e-01 + result * t;
              result = 1.7793648496568413e-03 + result * t;
              result = -3.8607931234291248e-01 + result * t;
              result = -1.3772243036256784e-03 + result * t;
              result = 2.4998361950946230e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -8.6211193418618600e-03;
              result = 3.9672659951507663e-02 + result * t;
              result = -2.8032356209313776e-02 + result * t;
              result = -1.2018100097997468e-01 + result * t;
              result = 1.5211454765440757e-01 + result * t;
              result = 1.8156741977359056e-01 + result * t;
              result = -2.8565580576367944e-01 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 3.8888850953375299e-03;
              result = -2.0675175441525351e-02 + result * t;
              result = 2.8960097320633170e-02 + result * t;
              result = 3.3007940280906339e-02 + result * t;
              result = -1.1721899629364071e-01 + result * t;
              result = 5.0547887857343300e-02 + result * t;
              result = 9.0625016096269703e-02 + result * t;
              result = -6.9135654915323988e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -1.0935783900081983e-03;
              result = 6.5470202258373581e-03 + result * t;
              result = -1.3424368326430808e-02 + result * t;
              result = 3.7917735980054659e-03 + result * t;
              result = 2.7021207542622858e-02 + result * t;
              result = -4.1921530752601251e-02 + result * t;
              result = 2.0067193675035979e-02 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.6864527449888091e-04;
              result = -1.1080285042200293e-03 + result * t;
              result = 2.8926068384211769e-03 + result * t;
              result = -3.4000082968751468e-03 + result * t;
              result = 6.0977953679685335e-04 + result * t;
              result = 2.8892073963802350e-03 + result * t;
              result = -3.0399198174633716e-03 + result * t;
              result = 9.8771757246140135e-04 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -1.0437674595874616e-05;
              result = 7.2488417272137105e-05 + result * t;
              result = -2.1401342242250001e-04 + result * t;
              result = 3.4518293939112904e-04 + result * t;
              result = -3.2217074343172047e-04 + result * t;
              result = 1.6568781090774194e-04 + result * t;
              result = -3.6819513535053763e-05 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 8.2186414140745018e-08;
              result = -5.7530489898521505e-07 + result * t;
              result = 1.7259146969556452e-06 + result * t;
              result = -2.8765244949260755e-06 + result * t;
              result = 2.8765244949260755e-06 + result * t;
              result = -1.7259146969556452e-06 + result * t;
              result = 5.7530489898521505e-07 + result * t;
              result = -8.2186414140745018e-08 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalNakSplineBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalNakSplineBasis: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return weaklyFundamentalSplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  WeaklyFundamentalSplineBasis<LT, IT> weaklyFundamentalSplineBasis;
};

// default type-def (unsigned int for level and index)
typedef WeaklyFundamentalNakSplineBasis<unsigned int, unsigned int> SWeaklyFundamentalNakSplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_HPP */
