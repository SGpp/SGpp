// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LAGRANGE_NOTAKNOT_SPLINE_BASE_DERIV3_HPP
#define LAGRANGE_NOTAKNOT_SPLINE_BASE_DERIV3_HPP

#include <sgpp/base/operation/hash/common/basis/LagrangeSplineBasisDeriv3.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Lagrange not-a-knot spline basis (3rd derivative).
 */
template <class LT, class IT>
class LagrangeNotAKnotSplineBasisDeriv3: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  LagrangeNotAKnotSplineBasisDeriv3() :
    lagrangeSplineBasisDeriv3(LagrangeSplineBasisDeriv3<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit LagrangeNotAKnotSplineBasisDeriv3(size_t degree) :
        lagrangeSplineBasisDeriv3(LagrangeSplineBasisDeriv3<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~LagrangeNotAKnotSplineBasisDeriv3() override {
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
    double innerDeriv = hInvDbl * hInvDbl * hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        return 0.0;

      case 3:
        if (l <= 1) {
          return 0.0;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return lagrangeSplineBasisDeriv3.eval(l, i, x);
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
              return innerDeriv * 1.2115384615384615e+00;
            } else {
              return innerDeriv * -3.4615384615384615e-01;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              return innerDeriv * 1.2857142857142858e+00;
            } else if (t < 2.0) {
              return innerDeriv * -7.5000000000000000e-01;
            } else {
              return innerDeriv * 1.0714285714285714e-01;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              return innerDeriv * 7.5000000000000000e-01;
            } else if (t < 0.0) {
              return innerDeriv * -3.5625000000000000e+00;
            } else if (t < 1.0) {
              return innerDeriv * 3.9375000000000000e+00;
            } else if (t < 2.0) {
              return innerDeriv * -1.7500000000000000e+00;
            } else {
              return innerDeriv * 2.5000000000000000e-01;
            }
          }
        }

      case 5:
        if (l <= 1) {
          return 0.0;
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return lagrangeSplineBasisDeriv3.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -4.0;
            result = 5.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 2.4687282282583428e-01;
              result = -1.2691424978379291e+00 + result * t;
              result = 1.5146032388775184e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.2372434806106664e-01;
              result = 2.1209443911707657e-01 + result * t;
              result = -7.0968849203760592e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.7750334637617993e-02;
              result = -3.5354257005056698e-02 + result * t;
              result = 1.7401241852249345e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.5837459443403190e-05;
              result = 1.4641227017928474e-04 + result * t;
              result = -2.0268051518936280e-04 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 8.3308627406534208e-01;
              result = -2.8730266276977825e+00 + result * t;
              result = 1.5203179094973753e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -2.0584475676024696e+00;
              result = 2.1254910166942702e+00 + result * t;
              result = 3.9901449299210745e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.1863623582426450e+00;
              result = -1.9914041185106688e+00 + result * t;
              result = 4.6605794208390827e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -8.7801456102014702e-02;
              result = 3.8132059797462137e-01 + result * t;
              result = -3.3898381818411549e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 2.4750562389288100e-01;
              result = -1.2705832064460740e+00 + result * t;
              result = 1.5139169463016338e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.2633382558747494e-01;
              result = 2.1445053691121194e-01 + result * t;
              result = -7.0282058000659287e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.9745509036264600e-02;
              result = -3.8217114263737938e-02 + result * t;
              result = 1.7834653323077702e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -6.3695190439563224e-04;
              result = 1.2739038087912645e-03 + result * t;
              result = -6.3695190439563224e-04 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 8.2796024981246197e-01;
              result = -2.8234797500053008e+00 + result * t;
              result = 1.4755483813198931e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -2.2588613385233618e+00;
              result = 2.1442817488694708e+00 + result * t;
              result = 4.5675137961614842e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.6930765637241338e+00;
              result = -2.3734409281772524e+00 + result * t;
              result = 3.4217178996225756e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -5.9333723333101451e-01;
              result = 1.0127121992710153e+00 + result * t;
              result = -3.3819257449086104e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 8.9880504818690407e-02;
              result = -1.7396226739101370e-01 + result * t;
              result = 8.1182391449139732e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.8993711231835614e-03;
              result = 5.7987422463671228e-03 + result * t;
              result = -2.8993711231835614e-03 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 2.1626503567880850e-01;
              result = -3.5781948689119553e-01 + result * t;
              result = -3.6362476740496488e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -3.1940746322270743e+00;
              result = 9.3977072718165544e-01 + result * t;
              result = 8.3656438369519337e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 7.0520428702544038e+00;
              result = -5.4483785372724931e+00 + result * t;
              result = -1.4177395213502255e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -7.7003182112619371e+00;
              result = 8.6557072032363145e+00 + result * t;
              result = 1.8592481163168528e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 4.7084312063360576e+00;
              result = -6.7449292192875596e+00 + result * t;
              result = 1.1413138036060626e+00 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.5639574866992620e+00;
              result = 2.6719331933845547e+00 + result * t;
              result = -8.9518420934544018e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 2.3559058634055086e-01;
              result = -4.5598178001396944e-01 + result * t;
              result = 2.1279149733985239e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -7.5996963335661569e-03;
              result = 1.5199392667132314e-02 + result * t;
              result = -7.5996963335661569e-03 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 7:
        if (l <= 1) {
          return 0.0;
        } else if ((i > 9) && (i < hInv - 9)) {
          // l >= 5, 9 < i < 2^l - 9
          return lagrangeSplineBasisDeriv3.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -4.0;
            result = 5.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 4.7257663953740496e-03;
              result = -7.2800774289728365e-02 + result * t;
              result = 4.0961853035033713e-01 + result * t;
              result = -9.8006357750919249e-01 + result * t;
              result = 8.1494671275333441e-01 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -2.4684069831375232e-04;
              result = 2.8114880362564193e-03 + result * t;
              result = -1.0317187170494553e-02 + result * t;
              result = 1.2243696602299558e-02 + result * t;
              result = -8.6446900490032101e-04 + result * t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 8.3133438877362933e-03;
              result = -1.0471170512353835e-01 + result * t;
              result = 4.2040688106992030e-01 + result * t;
              result = -5.2107233592401292e-01 + result * t;
              result = -1.4220019847024264e-02 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -3.8471994189810755e-03;
              result = 2.8301797080242358e-02 + result * t;
              result = -3.8052567189855620e-02 + result * t;
              result = -5.5763098033999758e-02 + result * t;
              result = 5.4667640929686186e-02 + result * t;
              return result;
            }
          } else if ((l == 4) && (i == 7)) {
            // l = 4, i = 7
            if ((t < -7.0) || (t > 9.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 9.3988868694092342e-03;
              result = -4.2569907009851832e-02 + result * t;
              result = 9.8017704228193815e-03 + result * t;
              result = 5.4297435637160225e-02 + result * t;
              result = 8.2879612814939236e-03 + result * t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -3.9757900880656227e-01;
              result = 1.0781228290069590e-01 + result * t;
              result = 4.0125602576788383e-01 + result * t;
              result = 4.9547110111559117e-01 + result * t;
              result = 6.3947020533491536e-02 + result * t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 1.3490940983998960e+00;
              result = -1.4825037523255533e+00 + result * t;
              result = -1.6607811783694022e+00 + result * t;
              result = 3.1103966127197402e-02 + result * t;
              result = 6.7090742151110017e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -2.3388884218038686e+00;
              result = 3.9138726412740308e+00 + result * t;
              result = 1.9862721550533140e+00 + result * t;
              result = -2.3415932539886830e+00 + result * t;
              result = -1.0921794446567619e+00 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 2.5086593678602669e+00;
              result = -5.4416810459414444e+00 + result * t;
              result = -3.0544045194780645e-01 + result * t;
              result = 4.0170152927245626e+00 + result * t;
              result = 1.2748367587803128e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.7464401069355677e+00;
              result = 4.5929564254996231e+00 + result * t;
              result = -1.5785273826105377e+00 + result * t;
              result = -2.8842712775543160e+00 + result * t;
              result = 9.0603683857360973e-01 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 7.7884079810446727e-01;
              result = -2.3928040022426473e+00 + result * t;
              result = 1.7217012522749267e+00 + result * t;
              result = 7.5178280598120839e-01 + result * t;
              result = -7.1024550302718803e-01 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.0588726900077955e-01;
              result = 7.2255919017522174e-01 + result * t;
              result = -7.8366596582621162e-01 + result * t;
              result = 1.3213649622098889e-01 + result * t;
              result = 1.4927535109076698e-01 + result * t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 2.5352087190622076e-02;
              result = -1.0098988582789646e-01 + result * t;
              result = 1.4868799069477640e-01 + result * t;
              result = -9.1066940908887187e-02 + result * t;
              result = 1.4417802659986517e-02 + result * t;
              return result;
            } else {
              t -= 5.0;
              double result = -2.9220989450600728e-05;
              result = 4.1846293459185186e-04 + result * t;
              result = -2.1691436451805155e-03 + result * t;
              result = 4.7477317594645354e-03 + result * t;
              result = -3.5989461913986573e-03 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 5.4677543571193584e-03;
              result = -7.8284460449090071e-02 + result * t;
              result = 4.0569386656526518e-01 + result * t;
              result = -8.8771986924832014e-01 + result * t;
              result = 6.7272228107034904e-01 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.7442661657090047e-03;
              result = 9.1996092648196671e-03 + result * t;
              result = -8.8152405403572193e-03 + result * t;
              result = -7.7922859966169595e-05 + result * t;
              result = 2.4843158021028195e-03 + result * t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 4.7361204782693723e-04;
              result = -1.7774553980163521e-03 + result * t;
              result = 2.3179902598477521e-03 + result * t;
              result = -1.0866408090576282e-03 + result * t;
              result = 4.6495500890091815e-05 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -2.9480326880967046e-05;
              result = 1.1699279329139679e-04 + result * t;
              result = -1.7270364723968098e-04 + result * t;
              result = 1.1142170789656837e-04 + result * t;
              result = -2.5998398509199284e-05 + result * t;
              return result;
            } else {
              t -= 6.0;
              double result = 2.3212855811785075e-07;
              result = -9.2851423247140300e-07 + result * t;
              result = 1.3927713487071044e-06 + result * t;
              result = -9.2851423247140300e-07 + result * t;
              result = 2.3212855811785075e-07 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 4.6426855925738086e-02;
              result = -5.1757067002077739e-01 + result * t;
              result = 1.8231421333866467e+00 + result * t;
              result = -1.9197939669328985e+00 + result * t;
              result = -2.1636992218147943e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.0534276349169057e-01;
              result = 2.2525902479103208e-01 + result * t;
              result = 6.9272262008175059e-02 + result * t;
              result = -2.9277394384808730e-01 + result * t;
              result = 3.5480579932473159e-02 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 6.7082448561148769e-02;
              result = -1.9611202917573020e-01 + result * t;
              result = 1.1299275543112787e-01 + result * t;
              result = 1.0017660057459678e-01 + result * t;
              result = -6.8104840608097572e-02 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.1216505803166939e-02;
              result = 7.2217765068864875e-02 + result * t;
              result = -7.2848640729170133e-02 + result * t;
              result = 6.1558181542569739e-03 + result * t;
              result = 1.6034934783045646e-02 + result * t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 3.3692730107623653e-03;
              result = -1.2648258143802884e-02 + result * t;
              result = 1.6505619658422850e-02 + result * t;
              result = -7.7541913101564341e-03 + result * t;
              result = 3.4337147383041997e-04 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -2.0885298651649839e-04;
              result = 8.2883389924657625e-04 + result * t;
              result = -1.2235167084116126e-03 + result * t;
              result = 7.8936561833007268e-04 + result * t;
              result = -1.8418531094368363e-04 + result * t;
              return result;
            } else {
              t -= 6.0;
              double result = 1.6445117048543181e-06;
              result = -6.5780468194172722e-06 + result * t;
              result = 9.8670702291259075e-06 + result * t;
              result = -6.5780468194172722e-06 + result * t;
              result = 1.6445117048543181e-06 + result * t;
              return result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 5.5054077563268536e-02;
              result = -4.3198608823597967e-01 + result * t;
              result = 7.9197489944718547e-01 + result * t;
              result = 2.1728035494955270e-01 + result * t;
              result = -5.5722396845922884e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -4.5688026973227647e-01;
              result = 4.4887915277631690e-01 + result * t;
              result = 8.9333328668920886e-01 + result * t;
              result = -8.8408828603242909e-02 + result * t;
              result = -5.6976994841200490e-01 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 7.1329790353061506e-01;
              result = -1.3786419261527889e+00 + result * t;
              result = -5.0131087337549929e-01 + result * t;
              result = 1.2173741241750196e+00 + result * t;
              result = 2.2715339271800147e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -5.8221141412178645e-01;
              result = 1.4745496879696711e+00 + result * t;
              result = -3.5744923065017614e-01 + result * t;
              result = -1.0679817869118857e+00 + result * t;
              result = 2.7787262089534792e-01 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 2.8384008058071114e-01;
              result = -8.5429596851747458e-01 + result * t;
              result = 5.7293134852811867e-01 + result * t;
              result = 3.1192315920962954e-01 + result * t;
              result = -2.5522012281882928e-01 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -8.2268812770407415e-02;
              result = 2.8106435380537004e-01 + result * t;
              result = -2.8691607354003817e-01 + result * t;
              result = 3.0258273036287725e-02 + result * t;
              result = 5.9178496982155479e-02 + result * t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.2788312030568892e-02;
              result = -4.8010897276259588e-02 + result * t;
              result = 6.2664111253627530e-02 + result * t;
              result = -2.9456063709308075e-02 + result * t;
              result = 1.3162375133676876e-03 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -7.9182253461117019e-04;
              result = 3.1423508460159824e-03 + result * t;
              result = -4.6387083917378792e-03 + result * t;
              result = 2.9927150914437926e-03 + result * t;
              result = -6.9830018800355164e-04 + result * t;
              return result;
            } else {
              t -= 6.0;
              double result = 6.2348231071745685e-06;
              result = -2.4939292428698274e-05 + result * t;
              result = 3.7408938643047408e-05 + result * t;
              result = -2.4939292428698274e-05 + result * t;
              result = 6.2348231071745685e-06 + result * t;
              return result;
            }
          } else if (i == 7) {
            // l >= 5, i = 7
            if ((t < -7.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 9.3943549117234768e-03;
              result = -4.2548117451891979e-02 + result * t;
              result = 9.7958660509280902e-03 + result * t;
              result = 5.4268697357702414e-02 + result * t;
              result = 8.2837753185900833e-03 + result * t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -3.9775314860395944e-01;
              result = 1.0776156113568365e-01 + result * t;
              result = 4.0107652815367806e-01 + result * t;
              result = 4.9528084547552204e-01 + result * t;
              result = 6.3967762044372350e-02 + result * t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 1.3515076848944523e+00;
              result = -1.4832510332801541e+00 + result * t;
              result = -1.6621576800630276e+00 + result * t;
              result = 2.9705990774091383e-02 + result * t;
              result = 6.7033354820529667e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -2.3495399524883211e+00;
              result = 3.9227797062976548e+00 + result * t;
              result = 1.9971353294632235e+00 + result * t;
              result = -2.3383317296146169e+00 + result * t;
              result = -1.0938614894693413e+00 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 2.5343114450509017e+00;
              result = -5.4753801036556293e+00 + result * t;
              result = -3.3176526657373800e-01 + result * t;
              result = 4.0261182382515104e+00 + result * t;
              result = 1.3818186418859887e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.7852866056941747e+00;
              result = 4.6618656765479773e+00 + result * t;
              result = -1.5520369072352156e+00 + result * t;
              result = -2.9263068256592462e+00 + result * t;
              result = 8.9146617726164379e-01 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 8.1760338853661418e-01;
              result = -2.4792807462287221e+00 + result * t;
              result = 1.7218404882436675e+00 + result * t;
              result = 8.1406996673755538e-01 + result * t;
              result = -7.1029848477901558e-01 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.3131980886368012e-01;
              result = 7.9113280791773466e-01 + result * t;
              result = -8.1038141922281359e-01 + result * t;
              result = 9.0322258685180984e-02 + result * t;
              result = 1.6393461251009947e-01 + result * t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 3.5730812040648631e-02;
              result = -1.3414642753698577e-01 + result * t;
              result = 1.7509815134830975e-01 + result * t;
              result = -8.2321391461962642e-02 + result * t;
              result = 3.6884510265214284e-03 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -2.2116194830403015e-03;
              result = 8.7768206256087552e-03 + result * t;
              result = -1.2956259018755783e-02 + result * t;
              result = 8.3588767862940536e-03 + result * t;
              result = -1.9504045834686124e-03 + result * t;
              return result;
            } else {
              t -= 6.0;
              double result = 1.7414326638112611e-05;
              result = -6.9657306552450444e-05 + result * t;
              result = 1.0448595982867566e-04 + result * t;
              result = -6.9657306552450444e-05 + result * t;
              result = 1.7414326638112611e-05 + result * t;
              return result;
            }
          } else {
            // l >= 5, i = 9
            if ((t < -9.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -5.0) {
              t += 9.0;
              double result = 2.8938603192514867e-05;
              result = -4.8598665186522662e-05 + result * t;
              result = -4.6740483585137805e-05 + result * t;
              result = 2.2793486061483184e-07 + result * t;
              result = 1.3134791083890904e-05 + result * t;
              return result;
            } else if (t < -4.0) {
              t += 5.0;
              double result = -2.5110444782711638e-02;
              result = 4.1441898589371519e-04 + result * t;
              result = 2.1481814406580173e-03 + result * t;
              result = 4.7018505545102307e-03 + result * t;
              result = 3.5641666385105007e-03 + result * t;
              return result;
            } else if (t < -3.0) {
              t += 4.0;
              double result = 2.0432375979851350e-01;
              result = -1.0002736014495284e-01 + result * t;
              result = -1.4727123029793066e-01 + result * t;
              result = -9.0200308737339141e-02 + result * t;
              result = -1.4281827163139174e-02 + result * t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -7.7819874007906809e-01;
              result = 7.1726767904910116e-01 + result * t;
              result = 7.7858924805829188e-01 + result * t;
              result = 1.3247018942599506e-01 + result * t;
              result = -1.4745696654484830e-01 + result * t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 1.7740967109416927e+00;
              result = -2.3955272812671713e+00 + result * t;
              result = -1.7388001552688133e+00 + result * t;
              result = 7.2865676237360999e-01 + result * t;
              result = 7.0267140990947174e-01 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -2.6413405384715198e+00;
              result = 4.7008595624996001e+00 + result * t;
              result = 1.7191982665798298e+00 + result * t;
              result = -2.8391385481987590e+00 + result * t;
              result = -9.2890255331120997e-01 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 2.6563054463918498e+00;
              result = -5.8645025913864783e+00 + result * t;
              result = -2.6266276750488347e-02 + result * t;
              result = 4.1364745185736220e+00 + result * t;
              result = 1.0676189097941048e-02 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.8104350617909906e+00;
              result = 4.7607191941809202e+00 + result * t;
              result = -1.6819413725588266e+00 + result * t;
              result = -2.8843440235193922e+00 + result * t;
              result = 9.1268728592644532e-01 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 8.1666587002088120e-01;
              result = -2.4810210529830421e+00 + result * t;
              result = 1.7376058392379901e+00 + result * t;
              result = 7.9219056674175203e-01 + result * t;
              result = -7.0331397776184423e-01 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.2965146190172162e-01;
              result = 7.8564242710048293e-01 + result * t;
              result = -8.0546209958584858e-01 + result * t;
              result = 9.1002566352131178e-02 + result * t;
              result = 1.6212724525573716e-01 + result * t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 3.5415507644764990e-02;
              result = -1.3296342050640353e-01 + result * t;
              result = 1.7355641030527061e-01 + result * t;
              result = -8.1600199125003520e-02 + result * t;
              result = 3.6586772207811203e-03 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -2.1919116651336694e-03;
              result = 8.6986100726564528e-03 + result * t;
              result = -1.2840805345350001e-02 + result * t;
              result = 8.2843905453870974e-03 + result * t;
              result = -1.9330244605903228e-03 + result * t;
              return result;
            } else {
              t -= 6.0;
              double result = 1.7259146969556454e-05;
              result = -6.9036587878225817e-05 + result * t;
              result = 1.0355488181733872e-04 + result * t;
              result = -6.9036587878225817e-05 + result * t;
              result = 1.7259146969556454e-05 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const {
    return lagrangeSplineBasisDeriv3.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  LagrangeSplineBasisDeriv3<LT, IT> lagrangeSplineBasisDeriv3;
};

// default type-def (unsigned int for level and index)
typedef LagrangeNotAKnotSplineBasisDeriv3<unsigned int, unsigned int>
SLagrangeNotAKnotSplineBaseDeriv3;

}  // namespace base
}  // namespace sgpp

#endif /* LAGRANGE_NOTAKNOT_SPLINE_BASE_DERIV3_HPP */
