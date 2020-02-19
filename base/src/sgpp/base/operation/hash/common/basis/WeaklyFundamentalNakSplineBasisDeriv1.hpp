// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_DERIV1_HPP
#define WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_DERIV1_HPP

#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv1.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental not-a-knot spline basis (1st derivative).
 */
template <class LT, class IT>
class WeaklyFundamentalNakSplineBasisDeriv1 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalNakSplineBasisDeriv1()
      : weaklyFundamentalSplineBasisDeriv1(WeaklyFundamentalSplineBasisDeriv1<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalNakSplineBasisDeriv1(size_t degree)
      : weaklyFundamentalSplineBasisDeriv1(WeaklyFundamentalSplineBasisDeriv1<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~WeaklyFundamentalNakSplineBasisDeriv1() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    double innerDeriv = hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        return weaklyFundamentalSplineBasisDeriv1.eval(l, i, x);

      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -1.0;
          } else {
            // l = 0, i = 1
            return 1.0;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * -2.0 * t;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return weaklyFundamentalSplineBasisDeriv1.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 6.0576923076923073e-01;
              result = -1.9038461538461537e+00 + result * t;
              result = 1.0961538461538463e+00 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -1.7307692307692307e-01;
              result = 5.1923076923076927e-01 + result * t;
              result = -2.8846153846153844e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 6.4285714285714290e-01;
              result = -1.9285714285714286e+00 + result * t;
              result = 1.0714285714285714e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -3.7500000000000000e-01;
              result = 6.4285714285714290e-01 + result * t;
              result = -2.1428571428571427e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = 5.3571428571428568e-02;
              result = -1.0714285714285714e-01 + result * t;
              result = 5.3571428571428568e-02 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 3.7500000000000000e-01;
              result = -3.7500000000000000e-01 + result * t;
              result = -1.2500000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.7812500000000000e+00;
              result = 1.1250000000000000e+00 + result * t;
              result = 6.2500000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 1.9687500000000000e+00;
              result = -2.4375000000000000e+00 + result * t;
              result = -3.1250000000000000e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -8.7500000000000000e-01;
              result = 1.5000000000000000e+00 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = 1.2500000000000000e-01;
              result = -2.5000000000000000e-01 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -1.0;
          } else {
            // l = 0, i = 1
            return 1.0;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * -2.0 * t;
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return weaklyFundamentalSplineBasisDeriv1.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -6.6666666666666663e-01;
            result = 2.5000000000000000e+00 + result * t;
            result = -1.6666666666666667e+00 + result * t;
            result = -8.3333333333333337e-01 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 2.0572735235486191e-02;
              result = -2.1152374963965487e-01 + result * t;
              result = 7.5730161943875918e-01 + result * t;
              result = -1.0538211192245104e+00 + result * t;
              result = 4.0130037316525197e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.0310362338422219e-02;
              result = 3.5349073186179428e-02 + result * t;
              result = -3.5484424601880296e-02 + result * t;
              result = 7.0276256987197698e-04 + result * t;
              result = 1.0801904244253702e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.4791945531348327e-03;
              result = -5.8923761675094494e-03 + result * t;
              result = 8.7006209261246727e-03 + result * t;
              result = -5.4603164290392076e-03 + result * t;
              result = 1.0589530600025915e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.1531216202835992e-06;
              result = 2.4402045029880789e-05 + result * t;
              result = -1.0134025759468140e-04 + result * t;
              result = 1.8057513322111785e-04 + result * t;
              result = -1.1392405728656110e-04 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 6.9423856172111845e-02;
              result = -4.7883777128296373e-01 + result * t;
              result = 7.6015895474868767e-01 + result * t;
              result = 3.4279594037898875e-01 + result * t;
              result = -6.2082201056206943e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -1.7153729730020581e-01;
              result = 3.5424850278237840e-01 + result * t;
              result = 1.9950724649605373e-01 + result * t;
              result = -5.2709368918082655e-01 + result * t;
              result = -5.6291071385875184e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 9.8863529853553761e-02;
              result = -3.3190068641844478e-01 + result * t;
              result = 2.3302897104195414e-01 + result * t;
              result = 2.4851712295759296e-01 + result * t;
              result = -2.0116630858847537e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -7.3167880085012252e-03;
              result = 6.3553432995770223e-02 + result * t;
              result = -1.6949190909205775e-01 + result * t;
              result = 1.1432712520038181e-01 + result * t;
              result = 4.7342628846180658e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 2.0625468657740082e-02;
              result = -2.1176386774101233e-01 + result * t;
              result = 7.5695847315081688e-01 + result * t;
              result = -1.0516045754277841e+00 + result * t;
              result = 3.9985284700395129e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.0527818798956244e-02;
              result = 3.5741756151868656e-02 + result * t;
              result = -3.5141029000329643e-02 + result * t;
              result = 7.2449505713247277e-05 + result * t;
              result = 1.0703911347564831e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.6454590863553833e-03;
              result = -6.3695190439563227e-03 + result * t;
              result = 8.9173266615388512e-03 + result * t;
              result = -5.0956152351650580e-03 + result * t;
              result = 8.4926920586084306e-04 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -5.3079325366302692e-05;
              result = 2.1231730146521077e-04 + result * t;
              result = -3.1847595219781612e-04 + result * t;
              result = 2.1231730146521077e-04 + result * t;
              result = -5.3079325366302692e-05 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 6.8996687484371821e-02;
              result = -4.7057995833421679e-01 + result * t;
              result = 7.3777419065994654e-01 + result * t;
              result = 3.3804341140746691e-01 + result * t;
              result = -6.0137181556895192e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -1.8823844487694680e-01;
              result = 3.5738029147824513e-01 + result * t;
              result = 2.2837568980807421e-01 + result * t;
              result = -4.8932807134454948e-01 + result * t;
              result = -6.4201054196767235e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.4108971364367781e-01;
              result = -3.9557348802954206e-01 + result * t;
              result = 1.7108589498112878e-01 + result * t;
              result = 2.8661040319854714e-01 + result * t;
              result = -1.5601158913194416e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -4.9444769444251205e-02;
              result = 1.6878536654516921e-01 + result * t;
              result = -1.6909628724543052e-01 + result * t;
              result = 6.4205836468897593e-03 + result * t;
              result = 4.7200934661867512e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 7.4900420682242012e-03;
              result = -2.8993711231835615e-02 + result * t;
              result = 4.0591195724569866e-02 + result * t;
              result = -2.3194968985468491e-02 + result * t;
              result = 3.8658281642447489e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.4161426026529680e-04;
              result = 9.6645704106118721e-04 + result * t;
              result = -1.4496855615917807e-03 + result * t;
              result = 9.6645704106118721e-04 + result * t;
              result = -2.4161426026529680e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 1.8022086306567375e-02;
              result = -5.9636581148532590e-02 + result * t;
              result = -1.8181238370248244e-02 + result * t;
              result = 5.6786695461489733e-02 + result * t;
              result = 2.9057441814890835e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -2.6617288601892286e-01;
              result = 1.5662845453027591e-01 + result * t;
              result = 4.1828219184759668e-01 + result * t;
              result = 2.8389689533889684e-01 + result * t;
              result = -1.1461231731129669e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 5.8767023918786698e-01;
              result = -9.0806308954541548e-01 + result * t;
              result = -7.0886976067511276e-01 + result * t;
              result = 5.2565509854922654e-01 + result * t;
              result = 4.7802233838654995e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -6.4169318427182809e-01;
              result = 1.4426178672060523e+00 + result * t;
              result = 9.2962405815842641e-02 + result * t;
              result = -1.2655927346857774e+00 + result * t;
              result = -2.5585174096884770e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 3.9236926719467141e-01;
              result = -1.1241548698812600e+00 + result * t;
              result = 5.7065690180303130e-01 + result * t;
              result = 6.8141294147675269e-01 + result * t;
              result = -3.9729082003259530e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.3032979055827182e-01;
              result = 4.4532219889742575e-01 + result * t;
              result = -4.4759210467272009e-01 + result * t;
              result = 1.9739204217720966e-02 + result * t;
              result = 1.2299342056060007e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.9632548861712572e-02;
              result = -7.5996963335661569e-02 + result * t;
              result = 1.0639574866992620e-01 + result * t;
              result = -6.0797570668529255e-02 + result * t;
              result = 1.0132928444754876e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -6.3330802779717978e-04;
              result = 2.5332321111887191e-03 + result * t;
              result = -3.7998481667830784e-03 + result * t;
              result = 2.5332321111887191e-03 + result * t;
              result = -6.3330802779717978e-04 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 7:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -1.0;
          } else {
            // l = 0, i = 1
            return 1.0;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * -2.0 * t;
        } else if ((i > 9) && (i < hInv - 9)) {
          // l >= 5, 9 < i < 2^l - 9
          return weaklyFundamentalSplineBasisDeriv1.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -6.6666666666666663e-01;
            result = 2.5000000000000000e+00 + result * t;
            result = -1.6666666666666667e+00 + result * t;
            result = -8.3333333333333337e-01 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.5752554651246830e-04;
              result = -3.6400387144864183e-03 + result * t;
              result = 3.4134877529194758e-02 + result * t;
              result = -1.6334392958486540e-01 + result * t;
              result = 4.0747335637666721e-01 + result * t;
              result = -4.7046751338266207e-01 + result * t;
              result = 1.6259926219922699e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -8.2280232771250766e-06;
              result = 1.4057440181282097e-04 + result * t;
              result = -8.5976559754121266e-04 + result * t;
              result = 2.0406161003832597e-03 + result * t;
              result = -4.3223450245016051e-04 + result * t;
              result = -4.0732317390158044e-03 + result * t;
              result = 2.6450596187038589e-03 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 2.7711146292454312e-04;
              result = -5.2355852561769174e-03 + result * t;
              result = 3.5033906755826690e-02 + result * t;
              result = -8.6845389320668820e-02 + result * t;
              result = -7.1100099235121318e-03 + result * t;
              result = 2.6864852561377756e-01 + result * t;
              result = -1.7219670451162619e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -1.2823998063270253e-04;
              result = 1.4150898540121178e-03 + result * t;
              result = -3.1710472658213014e-03 + result * t;
              result = -9.2938496723332931e-03 + result * t;
              result = 2.7333820464843093e-02 + result * t;
              result = 1.2893588627149366e-02 + result * t;
              result = -2.6978298050115802e-02 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 4) && (i == 7)) {
            // l = 4, i = 7
            if ((t < -7.0) || (t > 9.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 3.1329622898030777e-04;
              result = -2.1284953504925913e-03 + result * t;
              result = 8.1681420190161513e-04 + result * t;
              result = 9.0495726061933703e-03 + result * t;
              result = 4.1439806407469618e-03 + result * t;
              result = -1.1024417352188043e-02 + result * t;
              result = -6.7262882896154744e-03 + result * t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -1.3252633626885410e-02;
              result = 5.3906141450347951e-03 + result * t;
              result = 3.3438002147323655e-02 + result * t;
              result = 8.2578516852598524e-02 + result * t;
              result = 3.1973510266745768e-02 + result * t;
              result = -1.3397066921762324e-01 + result * t;
              result = -9.2561069964300122e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 4.4969803279996536e-02;
              result = -7.4125187616277663e-02 + result * t;
              result = -1.3839843153078352e-01 + result * t;
              result = 5.1839943545329004e-03 + result * t;
              result = 3.3545371075555008e-01 + result * t;
              result = 2.5890117942682001e-01 + result * t;
              result = -8.6403729397106027e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -7.7962947393462295e-02;
              result = 1.9569363206370155e-01 + result * t;
              result = 1.6552267958777617e-01 + result * t;
              result = -3.9026554233144711e-01 + result * t;
              result = -5.4608972232838093e-01 + result * t;
              result = 2.9095973947697573e-01 + result * t;
              result = 3.4558133927273232e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 8.3621978928675567e-02;
              result = -2.7208405229707222e-01 + result * t;
              result = -2.5453370995650537e-02 + result * t;
              result = 6.6950254878742710e-01 + result * t;
              result = 6.3741837939015639e-02 + result * t;
              result = -7.9923513786528899e-01 + result * t;
              result = -1.6560821652104624e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -5.8214670231185588e-02;
              result = 2.2964782127498118e-01 + result * t;
              result = -1.3154394855087814e-01 + result * t;
              result = -4.8071187959238598e-01 + result * t;
              result = 4.5301841928680486e-01 + result * t;
              result = 3.7625431247911362e-01 + result * t;
              result = -2.9646701715499812e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 2.5961359936815576e-02;
              result = -1.1964020011213236e-01 + result * t;
              result = 1.4347510435624389e-01 + result * t;
              result = 1.2529713433020140e-01 + result * t;
              result = -3.5512275151359401e-01 + result * t;
              result = 1.1293080305984521e-01 + result * t;
              result = 9.1983037511451843e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -6.8629089666926522e-03;
              result = 3.6127959508761090e-02 + result * t;
              result = -6.5305497152184297e-02 + result * t;
              result = 2.2022749370164817e-02 + result * t;
              result = 7.4637675545383489e-02 + result * t;
              result = -8.9955720491531443e-02 + result * t;
              result = 2.4884487568831534e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 8.4506957302073595e-04;
              result = -5.0494942913948231e-03 + result * t;
              result = 1.2390665891231365e-02 + result * t;
              result = -1.5177823484814531e-02 + result * t;
              result = 7.2089013299932586e-03 + result * t;
              result = 3.6282338446423307e-03 + result * t;
              result = -4.4512546172674602e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 5.0;
              double result = -9.7403298168669089e-07;
              result = 2.0923146729592594e-05 + result * t;
              result = -1.8076197043170964e-04 + result * t;
              result = 7.9128862657742256e-04 + result * t;
              result = -1.7994730956993287e-03 + result * t;
              result = 1.8981755962610195e-03 + result * t;
              result = -6.0570175458912344e-04 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.8225847857064527e-04;
              result = -3.9142230224545032e-03 + result * t;
              result = 3.3807822213772099e-02 + result * t;
              result = -1.4795331154138669e-01 + result * t;
              result = 3.3636114053517452e-01 + result * t;
              result = -3.5470301076655897e-01 + result * t;
              result = 1.1315257493880702e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -9.1475538856966820e-05;
              result = 4.5998046324098330e-04 + result * t;
              result = -7.3460337836310161e-04 + result * t;
              result = -1.2987143327694932e-05 + result * t;
              result = 1.2421579010514097e-03 + result * t;
              result = -1.1797301497864087e-03 + result * t;
              result = 2.7568174422425856e-04 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.5787068260897907e-05;
              result = -8.8872769900817604e-05 + result * t;
              result = 1.9316585498731266e-04 + result * t;
              result = -1.8110680150960470e-04 + result * t;
              result = 2.3247750445045908e-05 + result * t;
              result = 7.8259791944034938e-05 + result * t;
              result = -4.0976101817520529e-05 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -9.8267756269890155e-07;
              result = 5.8496396645698386e-06 + result * t;
              result = -1.4391970603306747e-05 + result * t;
              result = 1.8570284649428062e-05 + result * t;
              result = -1.2999199254599642e-05 + result * t;
              result = 4.4568683158627344e-06 + result * t;
              result = -4.9520759065141497e-07 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 7.7376186039283589e-09;
              result = -4.6425711623570150e-08 + result * t;
              result = 1.1606427905892537e-07 + result * t;
              result = -1.5475237207856716e-07 + result * t;
              result = 1.1606427905892537e-07 + result * t;
              result = -4.6425711623570150e-08 + result * t;
              result = 7.7376186039283589e-09 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 1.5475618641912696e-03;
              result = -2.5878533501038869e-02 + result * t;
              result = 1.5192851111555389e-01 + result * t;
              result = -3.1996566115548308e-01 + result * t;
              result = -1.0818496109073972e-01 + result * t;
              result = 1.0068791387240570e+00 + result * t;
              result = -5.8500272890037075e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -3.5114254497230191e-03;
              result = 1.1262951239551604e-02 + result * t;
              result = 5.7726885006812549e-03 + result * t;
              result = -4.8795657308014551e-02 + result * t;
              result = 1.7740289966236580e-02 + result * t;
              result = 6.0443772378156291e-02 + result * t;
              result = -3.3353929161460835e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 2.2360816187049592e-03;
              result = -9.8056014587865111e-03 + result * t;
              result = 9.4160629525939891e-03 + result * t;
              result = 1.6696100095766130e-02 + result * t;
              result = -3.4052420304048786e-02 + result * t;
              result = 7.8743378887307282e-03 + result * t;
              result = 9.5586901654273245e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -7.0721686010556466e-04;
              result = 3.6108882534432439e-03 + result * t;
              result = -6.0707200607641783e-03 + result * t;
              result = 1.0259696923761623e-03 + result * t;
              result = 8.0174673915228232e-03 + result * t;
              result = -8.0894682033952892e-03 + result * t;
              result = 1.9232509583878375e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.1230910035874551e-04;
              result = -6.3241290719014424e-04 + result * t;
              result = 1.3754683048685708e-03 + result * t;
              result = -1.2923652183594056e-03 + result * t;
              result = 1.7168573691520999e-04 + result * t;
              result = 5.5163552030496478e-04 + result * t;
              result = -2.8982882853496373e-04 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -6.9617662172166132e-06;
              result = 4.1441694962328812e-05 + result * t;
              result = -1.0195972570096772e-04 + result * t;
              result = 1.3156093638834544e-04 + result * t;
              result = -9.2092655471841814e-05 + result * t;
              result = 3.1574624733202907e-05 + result * t;
              result = -3.5082916370225450e-06 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 5.4817056828477266e-08;
              result = -3.2890234097086361e-07 + result * t;
              result = 8.2225585242715903e-07 + result * t;
              result = -1.0963411365695453e-06 + result * t;
              result = 8.2225585242715903e-07 + result * t;
              result = -3.2890234097086361e-07 + result * t;
              result = 5.4817056828477266e-08 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 1.8351359187756180e-03;
              result = -2.1599304411798984e-02 + result * t;
              result = 6.5997908287265455e-02 + result * t;
              result = 3.6213392491592118e-02 + result * t;
              result = -2.7861198422961442e-01 + result * t;
              result = -4.6392448883231088e-02 + result * t;
              result = 2.3267290767345314e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.5229342324409216e-02;
              result = 2.2443957638815844e-02 + result * t;
              result = 7.4444440557434072e-02 + result * t;
              result = -1.4734804767207151e-02 + result * t;
              result = -2.8488497420600245e-01 + result * t;
              result = -1.3615523729071554e-02 + result * t;
              result = 2.0146201109132153e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 2.3776596784353832e-02;
              result = -6.8932096307639454e-02 + result * t;
              result = -4.1775906114624936e-02 + result * t;
              result = 2.0289568736250327e-01 + result * t;
              result = 1.1357669635900074e-01 + result * t;
              result = -3.0896838996533771e-01 + result * t;
              result = -3.0114235739118922e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.9407047137392879e-02;
              result = 7.3727484398483548e-02 + result * t;
              result = -2.9787435887514677e-02 + result * t;
              result = -1.7799696448531430e-01 + result * t;
              result = 1.3893631044767396e-01 + result * t;
              result = 1.5776753954959963e-01 + result * t;
              result = -1.0954164762086316e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 9.4613360193570376e-03;
              result = -4.2714798425873730e-02 + result * t;
              result = 4.7744279044009887e-02 + result * t;
              result = 5.1987193201604928e-02 + result * t;
              result = -1.2761006140941464e-01 + result * t;
              result = 3.4694662607006439e-02 + result * t;
              result = 3.3698239264672129e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.7422937590135801e-03;
              result = 1.4053217690268502e-02 + result * t;
              result = -2.3909672795003178e-02 + result * t;
              result = 5.0430455060479539e-03 + result * t;
              result = 2.9589248491077740e-02 + result * t;
              result = -3.0392740444194943e-02 + result * t;
              result = 7.2608503013620491e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 4.2627706768562975e-04;
              result = -2.4005448638129796e-03 + result * t;
              result = 5.2220092711356275e-03 + result * t;
              result = -4.9093439515513456e-03 + result * t;
              result = 6.5811875668384380e-04 + result * t;
              result = 2.0885277733527097e-03 + result * t;
              result = -1.0983450094554580e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -2.6394084487039006e-05;
              result = 1.5711754230079912e-04 + result * t;
              result = -3.8655903264482323e-04 + result * t;
              result = 4.9878584857396551e-04 + result * t;
              result = -3.4915009400177582e-04 + result * t;
              result = 1.1970860365775171e-04 + result * t;
              result = -1.3300955961972412e-05 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 2.0782743690581893e-07;
              result = -1.2469646214349137e-06 + result * t;
              result = 3.1174115535872842e-06 + result * t;
              result = -4.1565487381163784e-06 + result * t;
              result = 3.1174115535872842e-06 + result * t;
              result = -1.2469646214349137e-06 + result * t;
              result = 2.0782743690581893e-07 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 7) {
            // l >= 5, i = 7
            if ((t < -7.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 3.1314516372411589e-04;
              result = -2.1274058725945992e-03 + result * t;
              result = 8.1632217091067414e-04 + result * t;
              result = 9.0447828929504023e-03 + result * t;
              result = 4.1418876592950417e-03 + result * t;
              result = -1.1018485378315195e-02 + result * t;
              result = -6.7227050766281628e-03 + result * t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -1.3258438286798648e-02;
              result = 5.3880780567841825e-03 + result * t;
              result = 3.3423044012806505e-02 + result * t;
              result = 8.2546807579253673e-02 + result * t;
              result = 3.1983881022186175e-02 + result * t;
              result = -1.3387096048932179e-01 + result * t;
              result = -9.2502886062100734e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 4.5050256163148408e-02;
              result = -7.4162551664007703e-02 + result * t;
              result = -1.3851314000525231e-01 + result * t;
              result = 4.9509984623485638e-03 + result * t;
              result = 3.3516677410264833e-01 + result * t;
              result = 2.5881916090716661e-01 + result * t;
              result = -8.6290474167190628e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -7.8317998416277360e-02;
              result = 1.9613898531488275e-01 + result * t;
              result = 1.6642794412193529e-01 + result * t;
              result = -3.8972195493576950e-01 + result * t;
              result = -5.4693074473467063e-01 + result * t;
              result = 2.8944192313735173e-01 + result * t;
              result = 3.4502102379886129e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 8.4477048168363392e-02;
              result = -2.7376900518278147e-01 + result * t;
              result = -2.7647105547811499e-02 + result * t;
              result = 6.7101970637525177e-01 + result * t;
              result = 6.9090932094299437e-02 + result * t;
              result = -7.9708671857480751e-01 + result * t;
              result = -1.7940821713686512e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -5.9509553523139164e-02;
              result = 2.3309328382739886e-01 + result * t;
              result = -1.2933640893626797e-01 + result * t;
              result = -4.8771780427654104e-01 + result * t;
              result = 4.4573308863082189e-01 + result * t;
              result = 3.8158310564457359e-01 + result * t;
              result = -2.9185596438117245e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 2.7253446284553805e-02;
              result = -1.2396403731143610e-01 + result * t;
              result = 1.4348670735363897e-01 + result * t;
              result = 1.3567832778959257e-01 + result * t;
              result = -3.5514924238950779e-01 + result * t;
              result = 1.0095933232968174e-01 + result * t;
              result = 9.1989746985673748e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -7.7106602954560037e-03;
              result = 3.9556640395886730e-02 + result * t;
              result = -6.7531784935234471e-02 + result * t;
              result = 1.5053709780863497e-02 + result * t;
              result = 8.1967306255049735e-02 + result * t;
              result = -8.4656848515857960e-02 + result * t;
              result = 2.0254281042196943e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.1910270680216210e-03;
              result = -6.7073213768492883e-03 + result * t;
              result = 1.4591512612359145e-02 + result * t;
              result = -1.3720231910327106e-02 + result * t;
              result = 1.8442255132607142e-03 + result * t;
              result = 5.8309938025917912e-03 + result * t;
              result = -3.0673562725515159e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -7.3720649434676723e-05;
              result = 4.3884103128043780e-04 + result * t;
              result = -1.0796882515629819e-03 + result * t;
              result = 1.3931461310490087e-03 + result * t;
              result = -9.7520229173430619e-04 + result * t;
              result = 3.3435507145176214e-04 + result * t;
              result = -3.7150563494640233e-05 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 5.8047755460375364e-07;
              result = -3.4828653276225220e-06 + result * t;
              result = 8.7071633190563055e-06 + result * t;
              result = -1.1609551092075074e-05 + result * t;
              result = 8.7071633190563055e-06 + result * t;
              result = -3.4828653276225220e-06 + result * t;
              result = 5.8047755460375364e-07 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 5, i = 9
            if ((t < -9.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -5.0) {
              t += 9.0;
              double result = 9.6462010641716217e-07;
              result = -2.4299332593261330e-06 + result * t;
              result = -3.8950402987614837e-06 + result * t;
              result = 3.7989143435805307e-08 + result * t;
              result = 6.5673955419454521e-06 + result * t;
              result = 6.2882884103033662e-06 + result * t;
              result = 1.4835950892407763e-06 + result * t;
              return innerDeriv * result;
            } else if (t < -4.0) {
              t += 5.0;
              double result = -8.3701482609038790e-04;
              result = 2.0720949294685758e-05 + result * t;
              result = 1.7901512005483479e-04 + result * t;
              result = 7.8364175908503838e-04 + result * t;
              result = 1.7820833192552503e-03 + result * t;
              result = 1.8798319770374401e-03 + result * t;
              result = 5.9984836443326931e-04 + result * t;
              return innerDeriv * result;
            } else if (t < -3.0) {
              t += 4.0;
              double result = 6.8107919932837834e-03;
              result = -5.0013680072476417e-03 + result * t;
              result = -1.2272602524827555e-02 + result * t;
              result = -1.5033384789556523e-02 + result * t;
              result = -7.1409135815695872e-03 + result * t;
              result = 3.5925001629534964e-03 + result * t;
              result = 4.4081266630701310e-03 + result * t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -2.5939958002635603e-02;
              result = 3.5863383952455058e-02 + result * t;
              result = 6.4882437338190990e-02 + result * t;
              result = 2.2078364904332508e-02 + result * t;
              result = -7.3728483272424150e-02 + result * t;
              result = -8.9021979544700980e-02 + result * t;
              result = -2.4636850083893897e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 5.9136557031389758e-02;
              result = -1.1977636406335856e-01 + result * t;
              result = -1.4490001293906776e-01 + result * t;
              result = 1.2144279372893499e-01 + result * t;
              result = 3.5133570495473587e-01 + result * t;
              result = 1.1296306972267388e-01 + result * t;
              result = -9.0503084708676074e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -8.8044684615717322e-02;
              result = 2.3504297812497998e-01 + result * t;
              result = 1.4326652221498581e-01 + result * t;
              result = -4.7318975803312652e-01 + result * t;
              result = -4.6445127665560498e-01 + result * t;
              result = 3.5630033093422525e-01 + result * t;
              result = 2.8969866372663206e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 8.8543514879728322e-02;
              result = -2.9322512956932395e-01 + result * t;
              result = -2.1888563958740287e-03 + result * t;
              result = 6.8941241976227030e-01 + result * t;
              result = 5.3380945489705242e-03 + result * t;
              result = -7.7215862468582497e-01 + result * t;
              result = -1.3772243036256784e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -6.0347835393033018e-02;
              result = 2.3803595970904598e-01 + result * t;
              result = -1.4016178104656887e-01 + result * t;
              result = -4.8072400391989872e-01 + result * t;
              result = 4.5634364296322266e-01 + result * t;
              result = 3.6313483954718112e-01 + result * t;
              result = -2.8565580576367944e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 2.7222195667362708e-02;
              result = -1.2405105264915210e-01 + result * t;
              result = 1.4480048660316586e-01 + result * t;
              result = 1.3203176112362536e-01 + result * t;
              result = -3.5165698888092212e-01 + result * t;
              result = 1.0109577571468660e-01 + result * t;
              result = 9.0625016096269703e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -7.6550487300573872e-03;
              result = 3.9282121355024145e-02 + result * t;
              result = -6.7121841632154039e-02 + result * t;
              result = 1.5167094392021864e-02 + result * t;
              result = 8.1063622627868578e-02 + result * t;
              result = -8.3843061505202501e-02 + result * t;
              result = 2.0067193675035979e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.1805169214921664e-03;
              result = -6.6481710253201761e-03 + result * t;
              result = 1.4463034192105883e-02 + result * t;
              result = -1.3600033187500587e-02 + result * t;
              result = 1.8293386103905602e-03 + result * t;
              result = 5.7784147927604701e-03 + result * t;
              result = -3.0399198174633716e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -7.3063722171122316e-05;
              result = 4.3493050363282261e-04 + result * t;
              result = -1.0700671121125001e-03 + result * t;
              result = 1.3807317575645162e-03 + result * t;
              result = -9.6651223029516140e-04 + result * t;
              result = 3.3137562181548389e-04 + result * t;
              result = -3.6819513535053763e-05 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 5.7530489898521505e-07;
              result = -3.4518293939112905e-06 + result * t;
              result = 8.6295734847782271e-06 + result * t;
              result = -1.1506097979704302e-05 + result * t;
              result = 8.6295734847782271e-06 + result * t;
              result = -3.4518293939112905e-06 + result * t;
              result = 5.7530489898521505e-07 + result * t;
              return innerDeriv * result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalNakSplineBasisDeriv1: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalNakSplineBasisDeriv1: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override {
    return weaklyFundamentalSplineBasisDeriv1.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  WeaklyFundamentalSplineBasisDeriv1<LT, IT> weaklyFundamentalSplineBasisDeriv1;
};

// default type-def (unsigned int for level and index)
typedef WeaklyFundamentalNakSplineBasisDeriv1<unsigned int, unsigned int>
    SWeaklyFundamentalNakSplineBaseDeriv1;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_DERIV1_HPP */
