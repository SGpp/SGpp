// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_DERIV2_HPP
#define WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_DERIV2_HPP

#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv2.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental not-a-knot spline basis (2nd derivative).
 */
template <class LT, class IT>
class WeaklyFundamentalNakSplineBasisDeriv2 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalNakSplineBasisDeriv2()
      : weaklyFundamentalSplineBasisDeriv2(WeaklyFundamentalSplineBasisDeriv2<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalNakSplineBasisDeriv2(size_t degree)
      : weaklyFundamentalSplineBasisDeriv2(WeaklyFundamentalSplineBasisDeriv2<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~WeaklyFundamentalNakSplineBasisDeriv2() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    double innerDeriv = hInvDbl * hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        return 0.0;

      case 3:
        if (l == 0) {
          return 0.0;
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * -2.0;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return weaklyFundamentalSplineBasisDeriv2.eval(l, i, x);
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
              double result = 1.2115384615384615e+00;
              result = -1.9038461538461537e+00 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -3.4615384615384615e-01;
              result = 5.1923076923076927e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.2857142857142858e+00;
              result = -1.9285714285714286e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -7.5000000000000000e-01;
              result = 6.4285714285714290e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = 1.0714285714285714e-01;
              result = -1.0714285714285714e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 7.5000000000000000e-01;
              result = -3.7500000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -3.5625000000000000e+00;
              result = 1.1250000000000000e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 3.9375000000000000e+00;
              result = -2.4375000000000000e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.7500000000000000e+00;
              result = 1.5000000000000000e+00 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = 2.5000000000000000e-01;
              result = -2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 0) {
          return 0.0;
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * -2.0;
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return weaklyFundamentalSplineBasisDeriv2.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -2.0;
            result = 5.0 + result * t;
            result = -1.6666666666666667e+00 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 8.2290940941944765e-02;
              result = -6.3457124891896455e-01 + result * t;
              result = 1.5146032388775184e+00 + result * t;
              result = -1.0538211192245104e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -4.1241449353688876e-02;
              result = 1.0604721955853828e-01 + result * t;
              result = -7.0968849203760592e-02 + result * t;
              result = 7.0276256987197698e-04 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 5.9167782125393309e-03;
              result = -1.7677128502528349e-02 + result * t;
              result = 1.7401241852249345e-02 + result * t;
              result = -5.4603164290392076e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -8.6124864811343967e-06;
              result = 7.3206135089642370e-05 + result * t;
              result = -2.0268051518936280e-04 + result * t;
              result = 1.8057513322111785e-04 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 2.7769542468844738e-01;
              result = -1.4365133138488912e+00 + result * t;
              result = 1.5203179094973753e+00 + result * t;
              result = 3.4279594037898875e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -6.8614918920082324e-01;
              result = 1.0627455083471351e+00 + result * t;
              result = 3.9901449299210745e-01 + result * t;
              result = -5.2709368918082655e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 3.9545411941421504e-01;
              result = -9.9570205925533439e-01 + result * t;
              result = 4.6605794208390827e-01 + result * t;
              result = 2.4851712295759296e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -2.9267152034004901e-02;
              result = 1.9066029898731068e-01 + result * t;
              result = -3.3898381818411549e-01 + result * t;
              result = 1.1432712520038181e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 8.2501874630960328e-02;
              result = -6.3529160322303702e-01 + result * t;
              result = 1.5139169463016338e+00 + result * t;
              result = -1.0516045754277841e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -4.2111275195824978e-02;
              result = 1.0722526845560597e-01 + result * t;
              result = -7.0282058000659287e-02 + result * t;
              result = 7.2449505713247277e-05 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 6.5818363454215332e-03;
              result = -1.9108557131868969e-02 + result * t;
              result = 1.7834653323077702e-02 + result * t;
              result = -5.0956152351650580e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.1231730146521077e-04;
              result = 6.3695190439563224e-04 + result * t;
              result = -6.3695190439563224e-04 + result * t;
              result = 2.1231730146521077e-04 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 2.7598674993748729e-01;
              result = -1.4117398750026504e+00 + result * t;
              result = 1.4755483813198931e+00 + result * t;
              result = 3.3804341140746691e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -7.5295377950778719e-01;
              result = 1.0721408744347354e+00 + result * t;
              result = 4.5675137961614842e-01 + result * t;
              result = -4.8932807134454948e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 5.6435885457471124e-01;
              result = -1.1867204640886262e+00 + result * t;
              result = 3.4217178996225756e-01 + result * t;
              result = 2.8661040319854714e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.9777907777700482e-01;
              result = 5.0635609963550765e-01 + result * t;
              result = -3.3819257449086104e-01 + result * t;
              result = 6.4205836468897593e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 2.9960168272896805e-02;
              result = -8.6981133695506849e-02 + result * t;
              result = 8.1182391449139732e-02 + result * t;
              result = -2.3194968985468491e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -9.6645704106118721e-04;
              result = 2.8993711231835614e-03 + result * t;
              result = -2.8993711231835614e-03 + result * t;
              result = 9.6645704106118721e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 7.2088345226269501e-02;
              result = -1.7890974344559776e-01 + result * t;
              result = -3.6362476740496488e-02 + result * t;
              result = 5.6786695461489733e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -1.0646915440756914e+00;
              result = 4.6988536359082772e-01 + result * t;
              result = 8.3656438369519337e-01 + result * t;
              result = 2.8389689533889684e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 2.3506809567514679e+00;
              result = -2.7241892686362466e+00 + result * t;
              result = -1.4177395213502255e+00 + result * t;
              result = 5.2565509854922654e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -2.5667727370873124e+00;
              result = 4.3278536016181572e+00 + result * t;
              result = 1.8592481163168528e-01 + result * t;
              result = -1.2655927346857774e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.5694770687786856e+00;
              result = -3.3724646096437798e+00 + result * t;
              result = 1.1413138036060626e+00 + result * t;
              result = 6.8141294147675269e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -5.2131916223308727e-01;
              result = 1.3359665966922774e+00 + result * t;
              result = -8.9518420934544018e-01 + result * t;
              result = 1.9739204217720966e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 7.8530195446850290e-02;
              result = -2.2799089000698472e-01 + result * t;
              result = 2.1279149733985239e-01 + result * t;
              result = -6.0797570668529255e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.5332321111887191e-03;
              result = 7.5996963335661569e-03 + result * t;
              result = -7.5996963335661569e-03 + result * t;
              result = 2.5332321111887191e-03 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 7:
        if (l == 0) {
          return 0.0;
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * -2.0;
        } else if ((i > 9) && (i < hInv - 9)) {
          // l >= 5, 9 < i < 2^l - 9
          return weaklyFundamentalSplineBasisDeriv2.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -2.0;
            result = 5.0 + result * t;
            result = -1.6666666666666667e+00 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 9.4515327907480982e-04;
              result = -1.8200193572432091e-02 + result * t;
              result = 1.3653951011677903e-01 + result * t;
              result = -4.9003178875459624e-01 + result * t;
              result = 8.1494671275333441e-01 + result * t;
              result = -4.7046751338266207e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -4.9368139662750459e-05;
              result = 7.0287200906410483e-04 + result * t;
              result = -3.4390623901648506e-03 + result * t;
              result = 6.1218483011497788e-03 + result * t;
              result = -8.6446900490032101e-04 + result * t;
              result = -4.0732317390158044e-03 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 1.6626687775472587e-03;
              result = -2.6177926280884587e-02 + result * t;
              result = 1.4013562702330676e-01 + result * t;
              result = -2.6053616796200646e-01 + result * t;
              result = -1.4220019847024264e-02 + result * t;
              result = 2.6864852561377756e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -7.6943988379621516e-04;
              result = 7.0754492700605895e-03 + result * t;
              result = -1.2684189063285205e-02 + result * t;
              result = -2.7881549016999879e-02 + result * t;
              result = 5.4667640929686186e-02 + result * t;
              result = 1.2893588627149366e-02 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 4) && (i == 7)) {
            // l = 4, i = 7
            if ((t < -7.0) || (t > 9.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 1.8797773738818467e-03;
              result = -1.0642476752462958e-02 + result * t;
              result = 3.2672568076064605e-03 + result * t;
              result = 2.7148717818580113e-02 + result * t;
              result = 8.2879612814939236e-03 + result * t;
              result = -1.1024417352188043e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -7.9515801761312455e-02;
              result = 2.6953070725173976e-02 + result * t;
              result = 1.3375200858929462e-01 + result * t;
              result = 2.4773555055779559e-01 + result * t;
              result = 6.3947020533491536e-02 + result * t;
              result = -1.3397066921762324e-01 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 2.6981881967997923e-01;
              result = -3.7062593808138833e-01 + result * t;
              result = -5.5359372612313407e-01 + result * t;
              result = 1.5551983063598701e-02 + result * t;
              result = 6.7090742151110017e-01 + result * t;
              result = 2.5890117942682001e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -4.6777768436077377e-01;
              result = 9.7846816031850770e-01 + result * t;
              result = 6.6209071835110467e-01 + result * t;
              result = -1.1707966269943415e+00 + result * t;
              result = -1.0921794446567619e+00 + result * t;
              result = 2.9095973947697573e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 5.0173187357205340e-01;
              result = -1.3604202614853611e+00 + result * t;
              result = -1.0181348398260215e-01 + result * t;
              result = 2.0085076463622813e+00 + result * t;
              result = 1.2748367587803128e-01 + result * t;
              result = -7.9923513786528899e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -3.4928802138711351e-01;
              result = 1.1482391063749058e+00 + result * t;
              result = -5.2617579420351257e-01 + result * t;
              result = -1.4421356387771580e+00 + result * t;
              result = 9.0603683857360973e-01 + result * t;
              result = 3.7625431247911362e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.5576815962089346e-01;
              result = -5.9820100056066183e-01 + result * t;
              result = 5.7390041742497555e-01 + result * t;
              result = 3.7589140299060420e-01 + result * t;
              result = -7.1024550302718803e-01 + result * t;
              result = 1.1293080305984521e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -4.1177453800155910e-02;
              result = 1.8063979754380544e-01 + result * t;
              result = -2.6122198860873719e-01 + result * t;
              result = 6.6068248110494443e-02 + result * t;
              result = 1.4927535109076698e-01 + result * t;
              result = -8.9955720491531443e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 5.0704174381244150e-03;
              result = -2.5247471456974115e-02 + result * t;
              result = 4.9562663564925460e-02 + result * t;
              result = -4.5533470454443593e-02 + result * t;
              result = 1.4417802659986517e-02 + result * t;
              result = 3.6282338446423307e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 5.0;
              double result = -5.8441978901201453e-06;
              result = 1.0461573364796296e-04 + result * t;
              result = -7.2304788172683856e-04 + result * t;
              result = 2.3738658797322677e-03 + result * t;
              result = -3.5989461913986573e-03 + result * t;
              result = 1.8981755962610195e-03 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.0935508714238716e-03;
              result = -1.9571115112272518e-02 + result * t;
              result = 1.3523128885508839e-01 + result * t;
              result = -4.4385993462416007e-01 + result * t;
              result = 6.7272228107034904e-01 + result * t;
              result = -3.5470301076655897e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -5.4885323314180097e-04;
              result = 2.2999023162049168e-03 + result * t;
              result = -2.9384135134524064e-03 + result * t;
              result = -3.8961429983084798e-05 + result * t;
              result = 2.4843158021028195e-03 + result * t;
              result = -1.1797301497864087e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 9.4722409565387443e-05;
              result = -4.4436384950408803e-04 + result * t;
              result = 7.7266341994925065e-04 + result * t;
              result = -5.4332040452881411e-04 + result * t;
              result = 4.6495500890091815e-05 + result * t;
              result = 7.8259791944034938e-05 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -5.8960653761934089e-06;
              result = 2.9248198322849197e-05 + result * t;
              result = -5.7567882413226987e-05 + result * t;
              result = 5.5710853948284183e-05 + result * t;
              result = -2.5998398509199284e-05 + result * t;
              result = 4.4568683158627344e-06 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 4.6425711623570150e-08;
              result = -2.3212855811785075e-07 + result * t;
              result = 4.6425711623570150e-07 + result * t;
              result = -4.6425711623570150e-07 + result * t;
              result = 2.3212855811785075e-07 + result * t;
              result = -4.6425711623570150e-08 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 9.2853711851476183e-03;
              result = -1.2939266750519435e-01 + result * t;
              result = 6.0771404446221555e-01 + result * t;
              result = -9.5989698346644925e-01 + result * t;
              result = -2.1636992218147943e-01 + result * t;
              result = 1.0068791387240570e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -2.1068552698338114e-02;
              result = 5.6314756197758019e-02 + result * t;
              result = 2.3090754002725020e-02 + result * t;
              result = -1.4638697192404365e-01 + result * t;
              result = 3.5480579932473159e-02 + result * t;
              result = 6.0443772378156291e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.3416489712229755e-02;
              result = -4.9028007293932550e-02 + result * t;
              result = 3.7664251810375957e-02 + result * t;
              result = 5.0088300287298389e-02 + result * t;
              result = -6.8104840608097572e-02 + result * t;
              result = 7.8743378887307282e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -4.2433011606333882e-03;
              result = 1.8054441267216219e-02 + result * t;
              result = -2.4282880243056713e-02 + result * t;
              result = 3.0779090771284870e-03 + result * t;
              result = 1.6034934783045646e-02 + result * t;
              result = -8.0894682033952892e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 6.7385460215247306e-04;
              result = -3.1620645359507209e-03 + result * t;
              result = 5.5018732194742834e-03 + result * t;
              result = -3.8770956550782170e-03 + result * t;
              result = 3.4337147383041997e-04 + result * t;
              result = 5.5163552030496478e-04 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -4.1770597303299677e-05;
              result = 2.0720847481164406e-04 + result * t;
              result = -4.0783890280387088e-04 + result * t;
              result = 3.9468280916503634e-04 + result * t;
              result = -1.8418531094368363e-04 + result * t;
              result = 3.1574624733202907e-05 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 3.2890234097086361e-07;
              result = -1.6445117048543181e-06 + result * t;
              result = 3.2890234097086361e-06 + result * t;
              result = -3.2890234097086361e-06 + result * t;
              result = 1.6445117048543181e-06 + result * t;
              result = -3.2890234097086361e-07 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 1.1010815512653708e-02;
              result = -1.0799652205899492e-01 + result * t;
              result = 2.6399163314906182e-01 + result * t;
              result = 1.0864017747477635e-01 + result * t;
              result = -5.5722396845922884e-01 + result * t;
              result = -4.6392448883231088e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -9.1376053946455288e-02;
              result = 1.1221978819407923e-01 + result * t;
              result = 2.9777776222973629e-01 + result * t;
              result = -4.4204414301621454e-02 + result * t;
              result = -5.6976994841200490e-01 + result * t;
              result = -1.3615523729071554e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 1.4265958070612300e-01;
              result = -3.4466048153819723e-01 + result * t;
              result = -1.6710362445849974e-01 + result * t;
              result = 6.0868706208750978e-01 + result * t;
              result = 2.2715339271800147e-01 + result * t;
              result = -3.0896838996533771e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.1644228282435728e-01;
              result = 3.6863742199241778e-01 + result * t;
              result = -1.1914974355005871e-01 + result * t;
              result = -5.3399089345594286e-01 + result * t;
              result = 2.7787262089534792e-01 + result * t;
              result = 1.5776753954959963e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 5.6768016116142232e-02;
              result = -2.1357399212936864e-01 + result * t;
              result = 1.9097711617603955e-01 + result * t;
              result = 1.5596157960481477e-01 + result * t;
              result = -2.5522012281882928e-01 + result * t;
              result = 3.4694662607006439e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -1.6453762554081482e-02;
              result = 7.0266088451342509e-02 + result * t;
              result = -9.5638691180012711e-02 + result * t;
              result = 1.5129136518143863e-02 + result * t;
              result = 5.9178496982155479e-02 + result * t;
              result = -3.0392740444194943e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 2.5576624061137785e-03;
              result = -1.2002724319064897e-02 + result * t;
              result = 2.0888037084542510e-02 + result * t;
              result = -1.4728031854654038e-02 + result * t;
              result = 1.3162375133676876e-03 + result * t;
              result = 2.0885277733527097e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -1.5836450692223404e-04;
              result = 7.8558771150399561e-04 + result * t;
              result = -1.5462361305792929e-03 + result * t;
              result = 1.4963575457218963e-03 + result * t;
              result = -6.9830018800355164e-04 + result * t;
              result = 1.1970860365775171e-04 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 1.2469646214349137e-06;
              result = -6.2348231071745685e-06 + result * t;
              result = 1.2469646214349137e-05 + result * t;
              result = -1.2469646214349137e-05 + result * t;
              result = 6.2348231071745685e-06 + result * t;
              result = -1.2469646214349137e-06 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 7) {
            // l >= 5, i = 7
            if ((t < -7.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 1.8788709823446953e-03;
              result = -1.0637029362972995e-02 + result * t;
              result = 3.2652886836426966e-03 + result * t;
              result = 2.7134348678851207e-02 + result * t;
              result = 8.2837753185900833e-03 + result * t;
              result = -1.1018485378315195e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -7.9550629720791882e-02;
              result = 2.6940390283920913e-02 + result * t;
              result = 1.3369217605122602e-01 + result * t;
              result = 2.4764042273776102e-01 + result * t;
              result = 6.3967762044372350e-02 + result * t;
              result = -1.3387096048932179e-01 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 2.7030153697889042e-01;
              result = -3.7081275832003852e-01 + result * t;
              result = -5.5405256002100922e-01 + result * t;
              result = 1.4852995387045691e-02 + result * t;
              result = 6.7033354820529667e-01 + result * t;
              result = 2.5881916090716661e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -4.6990799049766419e-01;
              result = 9.8069492657441371e-01 + result * t;
              result = 6.6571177648774116e-01 + result * t;
              result = -1.1691658648073084e+00 + result * t;
              result = -1.0938614894693413e+00 + result * t;
              result = 2.8944192313735173e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 5.0686228901018027e-01;
              result = -1.3688450259139073e+00 + result * t;
              result = -1.1058842219124600e-01 + result * t;
              result = 2.0130591191257552e+00 + result * t;
              result = 1.3818186418859887e-01 + result * t;
              result = -7.9708671857480751e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -3.5705732113883498e-01;
              result = 1.1654664191369943e+00 + result * t;
              result = -5.1734563574507186e-01 + result * t;
              result = -1.4631534128296231e+00 + result * t;
              result = 8.9146617726164379e-01 + result * t;
              result = 3.8158310564457359e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.6352067770732284e-01;
              result = -6.1982018655718052e-01 + result * t;
              result = 5.7394682941455588e-01 + result * t;
              result = 4.0703498336877769e-01 + result * t;
              result = -7.1029848477901558e-01 + result * t;
              result = 1.0095933232968174e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -4.6263961772736020e-02;
              result = 1.9778320197943366e-01 + result * t;
              result = -2.7012713974093788e-01 + result * t;
              result = 4.5161129342590492e-02 + result * t;
              result = 1.6393461251009947e-01 + result * t;
              result = -8.4656848515857960e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 7.1461624081297264e-03;
              result = -3.3536606884246443e-02 + result * t;
              result = 5.8366050449436581e-02 + result * t;
              result = -4.1160695730981321e-02 + result * t;
              result = 3.6884510265214284e-03 + result * t;
              result = 5.8309938025917912e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -4.4232389660806031e-04;
              result = 2.1942051564021888e-03 + result * t;
              result = -4.3187530062519276e-03 + result * t;
              result = 4.1794383931470268e-03 + result * t;
              result = -1.9504045834686124e-03 + result * t;
              result = 3.3435507145176214e-04 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 3.4828653276225220e-06;
              result = -1.7414326638112611e-05 + result * t;
              result = 3.4828653276225222e-05 + result * t;
              result = -3.4828653276225222e-05 + result * t;
              result = 1.7414326638112611e-05 + result * t;
              result = -3.4828653276225220e-06 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 5, i = 9
            if ((t < -9.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -5.0) {
              t += 9.0;
              double result = 5.7877206385029734e-06;
              result = -1.2149666296630666e-05 + result * t;
              result = -1.5580161195045935e-05 + result * t;
              result = 1.1396743030741592e-07 + result * t;
              result = 1.3134791083890904e-05 + result * t;
              result = 6.2882884103033662e-06 + result * t;
              return innerDeriv * result;
            } else if (t < -4.0) {
              t += 5.0;
              double result = -5.0220889565423281e-03;
              result = 1.0360474647342880e-04 + result * t;
              result = 7.1606048021933915e-04 + result * t;
              result = 2.3509252772551154e-03 + result * t;
              result = 3.5641666385105007e-03 + result * t;
              result = 1.8798319770374401e-03 + result * t;
              return innerDeriv * result;
            } else if (t < -3.0) {
              t += 4.0;
              double result = 4.0864751959702704e-02;
              result = -2.5006840036238209e-02 + result * t;
              result = -4.9090410099310221e-02 + result * t;
              result = -4.5100154368669571e-02 + result * t;
              result = -1.4281827163139174e-02 + result * t;
              result = 3.5925001629534964e-03 + result * t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -1.5563974801581362e-01;
              result = 1.7931691976227529e-01 + result * t;
              result = 2.5952974935276396e-01 + result * t;
              result = 6.6235094712997530e-02 + result * t;
              result = -1.4745696654484830e-01 + result * t;
              result = -8.9021979544700980e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 3.5481934218833855e-01;
              result = -5.9888182031679282e-01 + result * t;
              result = -5.7960005175627105e-01 + result * t;
              result = 3.6432838118680499e-01 + result * t;
              result = 7.0267140990947174e-01 + result * t;
              result = 1.1296306972267388e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -5.2826810769430388e-01;
              result = 1.1752148906249000e+00 + result * t;
              result = 5.7306608885994326e-01 + result * t;
              result = -1.4195692740993795e+00 + result * t;
              result = -9.2890255331120997e-01 + result * t;
              result = 3.5630033093422525e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 5.3126108927836990e-01;
              result = -1.4661256478466196e+00 + result * t;
              result = -8.7554255834961150e-03 + result * t;
              result = 2.0682372592868110e+00 + result * t;
              result = 1.0676189097941048e-02 + result * t;
              result = -7.7215862468582497e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -3.6208701235819807e-01;
              result = 1.1901797985452300e+00 + result * t;
              result = -5.6064712418627549e-01 + result * t;
              result = -1.4421720117596961e+00 + result * t;
              result = 9.1268728592644532e-01 + result * t;
              result = 3.6313483954718112e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.6333317400417624e-01;
              result = -6.2025526324576052e-01 + result * t;
              result = 5.7920194641266343e-01 + result * t;
              result = 3.9609528337087602e-01 + result * t;
              result = -7.0331397776184423e-01 + result * t;
              result = 1.0109577571468660e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -4.5930292380344323e-02;
              result = 1.9641060677512073e-01 + result * t;
              result = -2.6848736652861616e-01 + result * t;
              result = 4.5501283176065589e-02 + result * t;
              result = 1.6212724525573716e-01 + result * t;
              result = -8.3843061505202501e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 7.0831015289529991e-03;
              result = -3.3240855126600882e-02 + result * t;
              result = 5.7852136768423532e-02 + result * t;
              result = -4.0800099562501760e-02 + result * t;
              result = 3.6586772207811203e-03 + result * t;
              result = 5.7784147927604701e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -4.3838233302673392e-04;
              result = 2.1746525181641132e-03 + result * t;
              result = -4.2802684484500005e-03 + result * t;
              result = 4.1421952726935487e-03 + result * t;
              result = -1.9330244605903228e-03 + result * t;
              result = 3.3137562181548389e-04 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 3.4518293939112905e-06;
              result = -1.7259146969556454e-05 + result * t;
              result = 3.4518293939112908e-05 + result * t;
              result = -3.4518293939112908e-05 + result * t;
              result = 1.7259146969556454e-05 + result * t;
              result = -3.4518293939112905e-06 + result * t;
              return innerDeriv * result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalNakSplineBasisDeriv2: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalNakSplineBasisDeriv2: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override {
    return weaklyFundamentalSplineBasisDeriv2.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  WeaklyFundamentalSplineBasisDeriv2<LT, IT> weaklyFundamentalSplineBasisDeriv2;
};

// default type-def (unsigned int for level and index)
typedef WeaklyFundamentalNakSplineBasisDeriv2<unsigned int, unsigned int>
    SWeaklyFundamentalNakSplineBaseDeriv2;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_NAK_SPLINE_BASE_DERIV2_HPP */
