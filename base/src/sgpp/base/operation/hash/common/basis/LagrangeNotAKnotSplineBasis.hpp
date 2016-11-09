// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LAGRANGE_NOTAKNOT_SPLINE_BASE_HPP
#define LAGRANGE_NOTAKNOT_SPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/LagrangeSplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Lagrange not-a-knot spline basis.
 */
template <class LT, class IT>
class LagrangeNotAKnotSplineBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  LagrangeNotAKnotSplineBasis() : lagrangeSplineBasis(LagrangeSplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit LagrangeNotAKnotSplineBasis(size_t degree) :
        lagrangeSplineBasis(LagrangeSplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~LagrangeNotAKnotSplineBasis() override {
  }

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
        return lagrangeSplineBasis.eval(l, i, x);

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
          return lagrangeSplineBasis.eval(l, i, x);
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
              double result = 21.0/104.0;
              result = -99.0/104.0 + result * t;
              result = 57.0/52.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 1.0;
              double result = -3.0/52.0;
              result = 27.0/104.0 + result * t;
              result = -15.0/52.0 + result * t;
              result *= t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 3.0/14.0;
              result = -27.0/28.0 + result * t;
              result = 15.0/14.0 + result * t;
              result *= t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.0/8.0;
              result = 9.0/28.0 + result * t;
              result = -3.0/14.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = 1.0/56.0;
              result = -3.0/56.0 + result * t;
              result = 3.0/56.0 + result * t;
              result = -1.0/56.0 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 1.0/8.0;
              result = -3.0/16.0 + result * t;
              result = -1.0/8.0 + result * t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -19.0/32.0;
              result = 9.0/16.0 + result * t;
              result = 5.0/8.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 21.0/32.0;
              result = -39.0/32.0 + result * t;
              result = -1.0/32.0 + result * t;
              result = 19.0/32.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -7.0/24.0;
              result = 3.0/4.0 + result * t;
              result = -1.0/2.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = 1.0/24.0;
              result = -1.0/8.0 + result * t;
              result = 1.0/8.0 + result * t;
              result = -1.0/24.0 + result * t;
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
          return lagrangeSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -1.0/6.0;
            result = 5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 12096429.0/2939917532.0;
              result = -155465595.0/2939917532.0 + result * t;
              result = 556601077.0/2204938149.0 + result * t;
              result = -387268398.0/734979383.0 + result * t;
              result = 884842502.0/2204938149.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -6062323.0/2939917532.0;
              result = 6495210.0/734979383.0 + result * t;
              result = -52160641.0/4409876298.0 + result * t;
              result = 258258.0/734979383.0 + result * t;
              result = 95270123.0/8819752596.0 + result * t;
              result = -4483624.0/734979383.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 434871.0/1469958766.0;
              result = -4330775.0/2939917532.0 + result * t;
              result = 6394777.0/2204938149.0 + result * t;
              result = -2006610.0/734979383.0 + result * t;
              result = 2334926.0/2204938149.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -633.0/1469958766.0;
              result = 17935.0/2939917532.0 + result * t;
              result = -74483.0/2204938149.0 + result * t;
              result = 132719.0/1469958766.0 + result * t;
              result = -502391.0/4409876298.0 + result * t;
              result = 152131.0/2939917532.0 + result * t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 9840973451.0/708760215408.0;
              result = -5302830655.0/44297513463.0 + result * t;
              result = 11224383844.0/44297513463.0 + result * t;
              result = 399605468.0/2331448077.0 + result * t;
              result = -27500871371.0/44297513463.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = -8105254119.0/236253405136.0;
              result = 62769311285.0/708760215408.0 + result * t;
              result = 23567133167.0/354380107704.0 + result * t;
              result = -93395759171.0/354380107704.0 + result * t;
              result = -39896871881.0/708760215408.0 + result * t;
              result = 47033524987.0/236253405136.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 737584597.0/37303169232.0;
              result = -14702375125.0/177190053852.0 + result * t;
              result = 1146955998.0/14765837821.0 + result * t;
              result = 5504345300.0/44297513463.0 + result * t;
              result = -8911167263.0/44297513463.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = -1037169649.0/708760215408.0;
              result = 11261036215.0/708760215408.0 + result * t;
              result = -6673840111.0/118126702568.0 + result * t;
              result = 6752543157.0/118126702568.0 + result * t;
              result = 11184857273.0/236253405136.0 + result * t;
              result = -14750218887.0/236253405136.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 13665.0/3312652.0;
              result = -175375.0/3312652.0 + result * t;
              result = 89555.0/354927.0 + result * t;
              result = -435450.0/828163.0 + result * t;
              result = 993430.0/2484489.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -6975.0/3312652.0;
              result = 7400.0/828163.0 + result * t;
              result = -8315.0/709854.0 + result * t;
              result = 30.0/828163.0 + result * t;
              result = 106375.0/9937956.0 + result * t;
              result = -4850.0/828163.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 6541.0/19875912.0;
              result = -5275.0/3312652.0 + result * t;
              result = 1055.0/354927.0 + result * t;
              result = -2110.0/828163.0 + result * t;
              result = 2110.0/2484489.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -211.0/19875912.0;
              result = 1055.0/19875912.0 + result * t;
              result = -1055.0/9937956.0 + result * t;
              result = 1055.0/9937956.0 + result * t;
              result = -1055.0/19875912.0 + result * t;
              result = 211.0/19875912.0 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 630005.0/45654728.0;
              result = -2685525.0/22827364.0 + result * t;
              result = 601480.0/2445789.0 + result * t;
              result = 964580.0/5706841.0 + result * t;
              result = -10295800.0/17120523.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = -1718795.0/45654728.0;
              result = 4079025.0/45654728.0 + result * t;
              result = 744745.0/9783156.0 + result * t;
              result = -5585035.0/22827364.0 + result * t;
              result = -8793245.0/136964184.0 + result * t;
              result = 8265445.0/45654728.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 2576565.0/91309456.0;
              result = -2257475.0/22827364.0 + result * t;
              result = 139480.0/2445789.0 + result * t;
              result = 817820.0/5706841.0 + result * t;
              result = -2671000.0/17120523.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -902955.0/91309456.0;
              result = 3852925.0/91309456.0 + result * t;
              result = -7720045.0/136964184.0 + result * t;
              result = 146565.0/45654728.0 + result * t;
              result = 12929675.0/273928368.0 + result * t;
              result = -2406295.0/91309456.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 58621.0/39132624.0;
              result = -47275.0/6522104.0 + result * t;
              result = 66185.0/4891578.0 + result * t;
              result = -9455.0/815263.0 + result * t;
              result = 9455.0/2445789.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -1891.0/39132624.0;
              result = 9455.0/39132624.0 + result * t;
              result = -9455.0/19566312.0 + result * t;
              result = 9455.0/19566312.0 + result * t;
              result = -9455.0/39132624.0 + result * t;
              result = 1891.0/39132624.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 158239.0/43901410.0;
              result = -261813.0/17560564.0 + result * t;
              result = -57013.0/9407445.0 + result * t;
              result = 623254.0/21950705.0 + result * t;
              result = 1913494.0/65852115.0 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -2337073.0/43901410.0;
              result = 687621.0/17560564.0 + result * t;
              result = 9181589.0/65852115.0 + result * t;
              result = 6231737.0/43901410.0 + result * t;
              result = -15094927.0/131704230.0 + result * t;
              result = -1915181.0/12543260.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 134157671.0/1141436660.0;
              result = -3986525.0/17560564.0 + result * t;
              result = -15560191.0/65852115.0 + result * t;
              result = 1153850.0/4390141.0 + result * t;
              result = 31478782.0/65852115.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = -29298085.0/228287332.0;
              result = 41166423.0/114143666.0 + result * t;
              result = 53055349.0/1712154990.0 + result * t;
              result = -180574243.0/285359165.0 + result * t;
              result = -87611567.0/3424309980.0 + result * t;
              result = 45095679.0/114143666.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 537437599.0/6848619960.0;
              result = -64157579.0/228287332.0 + result * t;
              result = 162842177.0/856077495.0 + result * t;
              result = 13889102.0/40765595.0 + result * t;
              result = -9717478.0/24459357.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -25502263.0/978374280.0;
              result = 152492125.0/1369723992.0 + result * t;
              result = -39299849.0/263408460.0 + result * t;
              result = 33796577.0/3424309980.0 + result * t;
              result = 12959003.0/105363384.0 + result * t;
              result = -472077059.0/6848619960.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 22409311.0/5707183300.0;
              result = -2168643.0/114143666.0 + result * t;
              result = 1445762.0/40765595.0 + result * t;
              result = -8674572.0/285359165.0 + result * t;
              result = 2891524.0/285359165.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -722881.0/5707183300.0;
              result = 722881.0/1141436660.0 + result * t;
              result = -722881.0/570718330.0 + result * t;
              result = 722881.0/570718330.0 + result * t;
              result = -722881.0/1141436660.0 + result * t;
              result = 722881.0/5707183300.0 + result * t;
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
          return lagrangeSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -1.0/6.0;
            result = 5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 257939928953.0/11462137682716032.0;
              result = -3476885409655.0/5731068841358016.0 + result * t;
              result = 3260488883519.0/477589070113168.0 + result * t;
              result = -1125162626525.0/27553215583452.0 + result * t;
              result = 12162801335489.0/89547950646219.0 + result * t;
              result = -7021566944840.0/29849316882073.0 + result * t;
              result = 4853476902176.0/29849316882073.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 3.0;
              double result = -13472962237.0/11462137682716032.0;
              result = 44757865229.0/1910356280452672.0 + result * t;
              result = -82122930449.0/477589070113168.0 + result * t;
              result = 60910996615.0/119397267528292.0 + result * t;
              result = -4300634877.0/29849316882073.0 + result * t;
              result = -60791592456.0/29849316882073.0 + result * t;
              result = 236859668192.0/89547950646219.0 + result * t;
              result *= t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 74949333161255.0/1893264632909266656.0;
              result = -34417876382295.0/39443013185609722.0 + result * t;
              result = 1105474276890795.0/157772052742438888.0 + result * t;
              result = -1712721918042275.0/78886026371219444.0 + result * t;
              result = -280440215162905.0/118329039556829166.0 + result * t;
              result = 2649076834519710.0/19721506592804861.0 + result * t;
              result = -783687333065840.0/4551116906031891.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 1.0;
              double result = -34684602836695.0/1893264632909266656.0;
              result = 8586985810865.0/36408935248255128.0 + result * t;
              result = -100060527294385.0/157772052742438888.0 + result * t;
              result = -42297396365875.0/18204467624127564.0 + result * t;
              result = 359376080335965.0/39443013185609722.0 + result * t;
              result = 9780038196740.0/1517038968677297.0 + result * t;
              result = -1596158048574040.0/59164519778414583.0 + result * t;
              result *= t;
              return result;
            }
          } else if ((l == 4) && (i == 7)) {
            // l = 4, i = 7
            if ((t < -7.0) || (t > 9.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 62292654887911612529243.0/1391809233180362005872731520.0;
              result = -82290541155473880488063.0/231968205530060334312121920.0 + result * t;
              result = 56842477399975814540537.0/347952308295090501468182880.0 + result * t;
              result = 78720491935225032187315.0/34795230829509050146818288.0 + result * t;
              result = 360476907369518409416641.0/260964231221317876101137160.0 + result * t;
              result = -79916072193797486011057.0/14498012845628770894507620.0 + result * t;
              result = -219415081559056937471813.0/32620528902664734512642145.0 + result * t;
              result *= t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -276677067687834928775180669.0/146139969483938010616636809600.0;
              result = 625225544964340933945213.0/695904616590181002936365760.0 + result * t;
              result = 2326966006387491569426681.0/347952308295090501468182880.0 + result * t;
              result = 239444879620389088445243.0/11598410276503016715606096.0 + result * t;
              result = 2781314175402741163076401.0/260964231221317876101137160.0 + result * t;
              result = -323718080542430771692577.0/4832670948542923631502540.0 + result * t;
              result = -3019391058032024881386089.0/32620528902664734512642145.0 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 40819165708306763006082421.0/6353911716692956983332035200.0;
              result = -257920301338904700756824279.0/20877138497705430088090972800.0 + result * t;
              result = -192624214862224641349934269.0/6959046165901810029363657600.0 + result * t;
              result = 5411348405545321285776841.0/4175427699541086017618194560.0 + result * t;
              result = 466887571934189117599818947.0/4175427699541086017618194560.0 + result * t;
              result = 900852630018834225413613481.0/6959046165901810029363657600.0 + result * t;
              result = -1803862625341644565403914429.0/20877138497705430088090972800.0 + result * t;
              result = -17909619254706717312347840759.0/146139969483938010616636809600.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -12917803574782816224831413.0/1159841027650301671560609600.0;
              result = 1383984776325509854437137.0/42433208328669573349778400.0 + result * t;
              result = 1066555526625099566094041.0/32217806323619490876683600.0 + result * t;
              result = -1697422453778279816227513.0/17397615414754525073409144.0 + result * t;
              result = -293230215154915393226883.0/1610890316180974543834180.0 + result * t;
              result = 5272922550622493011077647.0/36245032114071927236269050.0 + result * t;
              result = 12525606739964169620035631.0/36245032114071927236269050.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 13855457424969562913863307.0/1159841027650301671560609600.0;
              result = -157787123411747332657614439.0/3479523082950905014681828800.0 + result * t;
              result = -5904372794551940214383729.0/1159841027650301671560609600.0 + result * t;
              result = 116477478630015843959131937.0/695904616590181002936365760.0 + result * t;
              result = 14786079763901377339818583.0/695904616590181002936365760.0 + result * t;
              result = -7855811047761926374829941.0/19658322502547485958654400.0 + result * t;
              result = -57623761211131186089489569.0/3479523082950905014681828800.0 + result * t;
              result = 308572798071966837499604077.0/1159841027650301671560609600.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -9645680420751674125845437.0/1159841027650301671560609600.0;
              result = 2080898164259585758336172.0/54367548171107890854403575.0 + result * t;
              result = -762850342342144608462079.0/28996025691257541789015240.0 + result * t;
              result = -261351262701628604580416.0/2174701926844315634176143.0 + result * t;
              result = 1641966715531634413292432.0/10873509634221578170880715.0 + result * t;
              result = 3409337409715881655229152.0/18122516057035963618134525.0 + result * t;
              result = -6447273934527611634185527.0/21747019268443156341761430.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 4301578626902221083312643.0/1159841027650301671560609600.0;
              result = -69381806323171668109239169.0/3479523082950905014681828800.0 + result * t;
              result = 11093887498585345361930893.0/386613675883433890520203200.0 + result * t;
              result = 21798713556476802476608247.0/695904616590181002936365760.0 + result * t;
              result = -27459062470501973317269979.0/231968205530060334312121920.0 + result * t;
              result = 65490889337152346170902031.0/1159841027650301671560609600.0 + result * t;
              result = 106685700753678557021823613.0/1159841027650301671560609600.0 + result * t;
              result = -85668439347383408533415093.0/1159841027650301671560609600.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -3494582951092364810250853.0/3564389499608244161381385600.0;
              result = 255504205387499690735687.0/42433208328669573349778400.0 + result * t;
              result = -92370725522199130384879.0/7072201388111595558296400.0 + result * t;
              result = 2336239779985700704777.0/424332083286695733497784.0 + result * t;
              result = 26392633629874202411567.0/1060830208216739333744460.0 + result * t;
              result = -39761560708049211707903.0/884025173513949444787050.0 + result * t;
              result = 65995540322526047047957.0/2652075520541848334361150.0 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 430308158930504769334939.0/3564389499608244161381385600.0;
              result = -428532486442368521422609.0/509198499944034880197340800.0 + result * t;
              result = 420620565674848638168539.0/169732833314678293399113600.0 + result * t;
              result = -386426247544145162330929.0/101839699988806976039468160.0 + result * t;
              result = 244717449565141677735899.0/101839699988806976039468160.0 + result * t;
              result = 307915205189675542618991.0/169732833314678293399113600.0 + result * t;
              result = -2266572173981550014548261.0/509198499944034880197340800.0 + result * t;
              result = 8095981205147039146749071.0/3564389499608244161381385600.0 + result * t;
              return result;
            } else {
              t -= 5.0;
              double result = -4723582220675034493.0/33946566662935658679822720.0;
              result = 59189082937874930411.0/16973283331467829339911360.0 + result * t;
              result = -306812413969181767087.0/8486641665733914669955680.0 + result * t;
              result = 2433109792729463435.0/12299480674976687927472.0 + result * t;
              result = -763574167516455736651.0/1272996249860087200493352.0 + result * t;
              result = 111869000722972813817.0/117870023135193259304940.0 + result * t;
              result = -481910038828518165761.0/795622656162554500308345.0 + result * t;
              result *= t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 17536589.0/673527640320.0;
              result = -219694783.0/336763820160.0 + result * t;
              result = 71157821.0/10523869380.0 + result * t;
              result = -103802755.0/2806365168.0 + result * t;
              result = 83289899.0/742861368.0 + result * t;
              result = -207380453.0/1169318820.0 + result * t;
              result = 893102189.0/7892902035.0 + result * t;
              result *= t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -52809689.0/4041165841920.0;
              result = 2868607.0/37418202240.0 + result * t;
              result = -773087.0/5261934690.0 + result * t;
              result = -27335.0/8419095504.0 + result * t;
              result = 5228923.0/12628643256.0 + result * t;
              result = -2069221.0/3507956460.0 + result * t;
              result = 2175929.0/7892902035.0 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 536119.0/237715637760.0;
              result = -59858267.0/4041165841920.0 + result * t;
              result = 52041017.0/1347055280640.0 + result * t;
              result = -36594131.0/808233168384.0 + result * t;
              result = 6263201.0/808233168384.0 + result * t;
              result = 52710133.0/1347055280640.0 + result * t;
              result = -165591223.0/4041165841920.0 + result * t;
              result = 53736667.0/4041165841920.0 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -189103.0/1347055280640.0;
              result = 72961.0/74836404480.0 + result * t;
              result = -323113.0/112254606720.0 + result * t;
              result = 52115.0/11225460672.0 + result * t;
              result = -72961.0/16838191008.0 + result * t;
              result = 10423.0/4677275280.0 + result * t;
              result = -10423.0/21047738760.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 1489.0/1347055280640.0;
              result = -10423.0/1347055280640.0 + result * t;
              result = 10423.0/449018426880.0 + result * t;
              result = -10423.0/269411056128.0 + result * t;
              result = 10423.0/269411056128.0 + result * t;
              result = -10423.0/449018426880.0 + result * t;
              result = 10423.0/1347055280640.0 + result * t;
              result = -1489.0/1347055280640.0 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 15006419435.0/67877697477312.0;
              result = -146381272345.0/33938848738656.0 + result * t;
              result = 64453484473.0/2121178046166.0 + result * t;
              result = -113117355995.0/1414118697444.0 + result * t;
              result = -229479564391.0/6363534138498.0 + result * t;
              result = 59326942339.0/117843224787.0 + result * t;
              result = -1861342418236.0/3181767069249.0 + result * t;
              result *= t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1021489175959.0/2036330924319360.0;
              result = 7078733305.0/3770983193184.0 + result * t;
              result = 2448980023.0/2121178046166.0 + result * t;
              result = -51752138515.0/4242356092332.0 + result * t;
              result = 37630313609.0/6363534138498.0 + result * t;
              result = 10684333583.0/353529674361.0 + result * t;
              result = -106124433436.0/3181767069249.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 650486021353.0/2036330924319360.0;
              result = -3327908247013.0/2036330924319360.0 + result * t;
              result = 1278281345047.0/678776974773120.0 + result * t;
              result = 1699939247027.0/407266184863872.0 + result * t;
              result = -4622799767537.0/407266184863872.0 + result * t;
              result = 116193896749.0/29512042381440.0 + result * t;
              result = 19464656379847.0/2036330924319360.0 + result * t;
              result = -14025154462459.0/2036330924319360.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -68577502973.0/678776974773120.0;
              result = 22694331527.0/37709831931840.0 + result * t;
              result = -68677749959.0/56564747897760.0 + result * t;
              result = 1450842925.0/5656474789776.0 + result * t;
              result = 22675301089.0/8484712184664.0 + result * t;
              result = -9532890199.0/2356864495740.0 + result * t;
              result = 20397788551.0/10605890230830.0 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 10890404483.0/678776974773120.0;
              result = -14308910665.0/135755394954624.0 + result * t;
              result = 210991235.0/766979632512.0 + result * t;
              result = -43861387661.0/135755394954624.0 + result * t;
              result = 7769088341.0/135755394954624.0 + result * t;
              result = 12481249655.0/45251798318208.0 + result * t;
              result = -39345827087.0/135755394954624.0 + result * t;
              result = 63898793077.0/678776974773120.0 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -18751931.0/18854915965920.0;
              result = 7234997.0/1047495331440.0 + result * t;
              result = -32040701.0/1571242997160.0 + result * t;
              result = 5167855.0/157124299716.0 + result * t;
              result = -7234997.0/235686449574.0 + result * t;
              result = 1033571.0/65468458215.0 + result * t;
              result = -2067142.0/589216123935.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 147653.0/18854915965920.0;
              result = -1033571.0/18854915965920.0 + result * t;
              result = 1033571.0/6284971988640.0 + result * t;
              result = -1033571.0/3770983193184.0 + result * t;
              result = 1033571.0/3770983193184.0 + result * t;
              result = -1033571.0/6284971988640.0 + result * t;
              result = 1033571.0/18854915965920.0 + result * t;
              result = -147653.0/18854915965920.0 + result * t;
              return result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 5163736075451513.0/19696716825354766080.0;
              result = -35452948551986707.0/9848358412677383040.0 + result * t;
              result = 4062319095625007.0/307761200396168220.0 + result * t;
              result = 19051413920735.0/2104350088178928.0 + result * t;
              result = -171491917422528829.0/1846567202377009320.0 + result * t;
              result = -264403625141131.0/11398562977635860.0 + result * t;
              result = 53705770023936569.0/230820900297126165.0 + result * t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1285577328001193413.0/590901504760642982400.0;
              result = 188919776945305.0/50504402116294272.0 + result * t;
              result = 229111103887771.0/15388060019808411.0 + result * t;
              result = -906960240551771.0/246208960316934576.0 + result * t;
              result = -175353083272941469.0/1846567202377009320.0 + result * t;
              result = -698388321146929.0/102587066798722740.0 + result * t;
              result = 3577049444289881.0/17755453869009705.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 2007089545423106491.0/590901504760642982400.0;
              result = -6788679905748285391.0/590901504760642982400.0 + result * t;
              result = -1645696385724748091.0/196967168253547660800.0 + result * t;
              result = 5994568348597407929.0/118180300952128596480.0 + result * t;
              result = 4474176052285078573.0/118180300952128596480.0 + result * t;
              result = -2340639571205015227.0/15151320634888281600.0 + result * t;
              result = -1368811324074008087.0/45453961904664844800.0 + result * t;
              result = 66454448003044689647.0/590901504760642982400.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -23742553533014869.0/8563789924067289600.0;
              result = 403385939567414447.0/32827861375591276800.0 + result * t;
              result = -2507327733460627.0/420870017635785600.0 + result * t;
              result = -14608149188499853.0/328278613755912768.0 + result * t;
              result = 38008182828436289.0/820696534389781920.0 + result * t;
              result = 53949697061483027.0/683913778658151600.0 + result * t;
              result = -112375563217236361.0/1025870667987227400.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 266224651946864093.0/196967168253547660800.0;
              result = -1402235481410907227.0/196967168253547660800.0 + result * t;
              result = 626937029573721643.0/65655722751182553600.0 + result * t;
              result = 511988511518510341.0/39393433650709532160.0 + result * t;
              result = -1675666162431581621.0/39393433650709532160.0 + result * t;
              result = 1138951574535717733.0/65655722751182553600.0 + result * t;
              result = 6637446763092982123.0/196967168253547660800.0 + result * t;
              result = -383131807030146209.0/15151320634888281600.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -494635381165093.0/1262610052907356800.0;
              result = 9611189212857113.0/4103482671948909600.0 + result * t;
              result = -490564640031319.0/102587066798722740.0 + result * t;
              result = 25867562309897.0/20517413359744548.0 + result * t;
              result = 1517737105739099.0/153880600198084110.0 + result * t;
              result = -649562936280181.0/42744611166134475.0 + result * t;
              result = 558652001161016.0/76940300099042055.0 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 999554606112671.0/16413930687795638400.0;
              result = -6567062834595011.0/16413930687795638400.0 + result * t;
              result = 745336506325621.0/713649160338940800.0 + result * t;
              result = -1343027189055209.0/1094262045853042560.0 + result * t;
              result = 2160463131309371.0/9848358412677383040.0 + result * t;
              result = 5713491685224589.0/5471310229265212800.0 + result * t;
              result = -54084476569464397.0/49241792063386915200.0 + result * t;
              result = 5857088334119939.0/16413930687795638400.0 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -185670288573463.0/49241792063386915200.0;
              result = 71636568032281.0/2735655114632606400.0 + result * t;
              result = -317247658428673.0/4103482671948909600.0 + result * t;
              result = 10233795433183.0/82069653438978192.0 + result * t;
              result = -71636568032281.0/615522400792336440.0 + result * t;
              result = 10233795433183.0/170978444664537900.0 + result * t;
              result = -10233795433183.0/769403000990420550.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 1461970776169.0/49241792063386915200.0;
              result = -10233795433183.0/49241792063386915200.0 + result * t;
              result = 10233795433183.0/16413930687795638400.0 + result * t;
              result = -10233795433183.0/9848358412677383040.0 + result * t;
              result = 10233795433183.0/9848358412677383040.0 + result * t;
              result = -10233795433183.0/16413930687795638400.0 + result * t;
              result = 10233795433183.0/49241792063386915200.0 + result * t;
              result = -1461970776169.0/49241792063386915200.0 + result * t;
              return result;
            }
          } else if (i == 7) {
            // l >= 5, i = 7
            if ((t < -7.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 358629707999.0/8016754677408000.0;
              result = -52638552407.0/148458419952000.0 + result * t;
              result = 40901591137.0/250523583669000.0 + result * t;
              result = 45318628477.0/20041886693520.0 + result * t;
              result = 1037640539561.0/751570751007000.0 + result * t;
              result = -46006507393.0/8350786122300.0 + result * t;
              result = -126314712581.0/18789268775175.0 + result * t;
              result *= t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -455527059217793.0/240502640322240000.0;
              result = 3599574996997.0/4008377338704000.0 + result * t;
              result = 1674652152643.0/250523583669000.0 + result * t;
              result = 689330735173.0/33403144489200.0 + result * t;
              result = 8012716493321.0/751570751007000.0 + result * t;
              result = -931606465861.0/13917976870500.0 + result * t;
              result = -8690307943501.0/93946343875875.0 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 1547815079204351.0/240502640322240000.0;
              result = -2972714914704731.0/240502640322240000.0 + result * t;
              result = -246761302893239.0/8907505197120000.0 + result * t;
              result = 59536410121309.0/48100528064448000.0 + result * t;
              result = 597099956444257.0/5344503118272000.0 + result * t;
              result = 10374448594026709.0/80167546774080000.0 + result * t;
              result = -56241427837039.0/651768672960000.0 + result * t;
              result = -29449983987193013.0/240502640322240000.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -896937400184177.0/80167546774080000.0;
              result = 1310331773287621.0/40083773387040000.0 + result * t;
              result = 667105999745461.0/20041886693520000.0 + result * t;
              result = -13017938771333.0/133612577956800.0 + result * t;
              result = -548076200758739.0/3006283004028000.0 + result * t;
              result = 40284459915787.0/278359537410000.0 + result * t;
              result = 1296538549848571.0/3757853755035000.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 967473958624783.0/80167546774080000.0;
              result = -1219299418237999.0/26722515591360000.0 + result * t;
              result = -443280125434421.0/80167546774080000.0 + result * t;
              result = 65602443533131.0/391061203776000.0 + result * t;
              result = 1107770106066907.0/48100528064448000.0 + result * t;
              result = -10650081132390637.0/26722515591360000.0 + result * t;
              result = -4314814991692181.0/240502640322240000.0 + result * t;
              result = 21227435048189573.0/80167546774080000.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -170383389841817.0/20041886693520000.0;
              result = 778604863914871.0/20041886693520000.0 + result * t;
              result = -51842913064949.0/2004188669352000.0 + result * t;
              result = -24436962429307.0/200418866935200.0 + result * t;
              result = 223333301447291.0/1503141502014000.0 + result * t;
              result = 79662972578021.0/417539306115000.0 + result * t;
              result = -535000991063.0/1833099392700.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 78030068920423.0/20041886693520000.0;
              result = -51759858122231.0/2505235836690000.0 + result * t;
              result = 575148866161579.0/20041886693520000.0 + result * t;
              result = 5665103484011.0/167015722446000.0 + result * t;
              result = -1423572175051997.0/12025132016112000.0 + result * t;
              result = 42154489566769.0/835078612230000.0 + result * t;
              result = 5530944258157339.0/60125660080560000.0 + result * t;
              result = -175206004010881.0/2505235836690000.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -231804269960629.0/210439810281960000.0;
              result = 4893763609819.0/742292099760000.0 + result * t;
              result = -135346438188313.0/10020943346760000.0 + result * t;
              result = 754261864363.0/200418866935200.0 + result * t;
              result = 41069486613419.0/1503141502014000.0 + result * t;
              result = -17673780893597.0/417539306115000.0 + result * t;
              result = 38056313034977.0/1878926877517500.0 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 35805644319307.0/210439810281960000.0;
              result = -67213687525919.0/60125660080560000.0 + result * t;
              result = 29244144246397.0/10020943346760000.0 + result * t;
              result = -41246900003339.0/12025132016112000.0 + result * t;
              result = 3696175877407.0/6012566008056000.0 + result * t;
              result = 58432058551081.0/20041886693520000.0 + result * t;
              result = -92213410294703.0/30062830040280000.0 + result * t;
              result = 419454168945661.0/420879620563920000.0 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -591000361177.0/56117282741856000.0;
              result = 32574823057.0/445375259856000.0 + result * t;
              result = -144259930681.0/668062889784000.0 + result * t;
              result = 4653546151.0/13361257795680.0 + result * t;
              result = -32574823057.0/100209433467600.0 + result * t;
              result = 4653546151.0/27835953741000.0 + result * t;
              result = -4653546151.0/125261791834500.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 4653546151.0/56117282741856000.0;
              result = -4653546151.0/8016754677408000.0 + result * t;
              result = 4653546151.0/2672251559136000.0 + result * t;
              result = -4653546151.0/1603350935481600.0 + result * t;
              result = 4653546151.0/1603350935481600.0 + result * t;
              result = -4653546151.0/2672251559136000.0 + result * t;
              result = 4653546151.0/8016754677408000.0 + result * t;
              result = -4653546151.0/56117282741856000.0 + result * t;
              return result;
            }
          } else {
            // l >= 5, i = 9
            if ((t < -9.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < -5.0) {
              t += 9.0;
              double result = 101069982853.0/733438869109615080.0;
              result = -84050533.0/207537880336620.0 + result * t;
              result = -142838697593.0/183359717277403770.0 + result * t;
              result = 174141965.0/18335971727740377.0 + result * t;
              result = 6337872578.0/2895153430695849.0 + result * t;
              result = 96084898756.0/30559952879567295.0 + result * t;
              result = 408047364176.0/275039575916105655.0 + result * t;
              result *= t;
              return result;
            } else if (t < -4.0) {
              t += 5.0;
              double result = -29467161958834333.0/246435460020830666880.0;
              result = 1266462468131.0/366719434554807540.0 + result * t;
              result = 211768785817.0/5914829589593670.0 + result * t;
              result = 1197402761605.0/6111990575913459.0 + result * t;
              result = 32676229358342.0/55007915183221131.0 + result * t;
              result = 9574596106628.0/10186650959855765.0 + result * t;
              result = 164982039767696.0/275039575916105655.0 + result * t;
              result *= t;
              return result;
            } else if (t < -3.0) {
              t += 4.0;
              double result = 239774379710154211.0/246435460020830666880.0;
              result = -29345581561893757.0/35205065717261523840.0 + result * t;
              result = -9601283964630751.0/3911673968584613760.0 + result * t;
              result = -26462564973460861.0/7041013143452304768.0 + result * t;
              result = -5586585153787615.0/2347004381150768256.0 + result * t;
              result = 21079034054341763.0/11735021905753841280.0 + result * t;
              result = 51729462954465569.0/11735021905753841280.0 + result * t;
              result = 554328359028715139.0/246435460020830666880.0 + result * t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -14495522637869689.0/3911673968584613760.0;
              result = 35071466358043409.0/5867510952876920640.0 + result * t;
              result = 2003675851216769.0/154408182970445280.0 + result * t;
              result = 539771032907435.0/97791849214615344.0 + result * t;
              result = -10815067078498789.0/440063321465769048.0 + result * t;
              result = -1813671666754577.0/40746603839423060.0 + result * t;
              result = -13552217597965699.0/550079151832211310.0 + result * t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 99138398884888597.0/11735021905753841280.0;
              result = -234263042679176651.0/11735021905753841280.0 + result * t;
              result = -68016193039359011.0/2347004381150768256.0 + result * t;
              result = 71256692235249869.0/2347004381150768256.0 + result * t;
              result = 824586438783458729.0/7041013143452304768.0 + result * t;
              result = 220937349622796069.0/3911673968584613760.0 + result * t;
              result = -637233408956765555.0/7041013143452304768.0 + result * t;
              result = -856144950764709733.0/11735021905753841280.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -14760090037866157.0/1173502190575384128.0;
              result = 57463218689380441.0/1466877738219230160.0 + result * t;
              result = 42030894413850697.0/1466877738219230160.0 + result * t;
              result = -17352788050303435.0/146687773821923016.0 + result * t;
              result = -34064661910680389.0/220031660732884524.0 + result * t;
              result = 10888521324324157.0/61119905759134590.0 + result * t;
              result = 79678597614835387.0/275039575916105655.0 + result * t;
              result *= t;
              return result;
            } else if (t < 1.0) {
              double result = 74218577623289503.0/5867510952876920640.0;
              result = -286750276567793731.0/5867510952876920640.0 + result * t;
              result = -2568627775413113.0/5867510952876920640.0 + result * t;
              result = 67418915400075035.0/391167396858461376.0 + result * t;
              result = 202073085377917.0/113564728120198464.0 + result * t;
              result = -755108197950409531.0/1955836984292306880.0 + result * t;
              result = -1275928213593467.0/926449097822671680.0 + result * t;
              result = 293356325102317331.0/1173502190575384128.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -151753536493300589.0/17602532858630761920.0;
              result = 250300824511003.0/6309151562233248.0 + result * t;
              result = -2570002454579624.0/91679858638701885.0 + result * t;
              result = -17629083489442835.0/146687773821923016.0 + result * t;
              result = 8367504135507700.0/55007915183221131.0 + result * t;
              result = 2219476717098217.0/12223981151826918.0 + result * t;
              result = -78566651675215844.0/275039575916105655.0 + result * t;
              result *= t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 3602854088111489.0/926449097822671680.0;
              result = -363935455067405753.0/17602532858630761920.0 + result * t;
              result = 169923688225196681.0/5867510952876920640.0 + result * t;
              result = 116204670678075139.0/3520506571726152384.0 + result * t;
              result = -412670246782905619.0/3520506571726152384.0 + result * t;
              result = 296590285647756127.0/5867510952876920640.0 + result * t;
              result = 1595229823648529201.0/17602532858630761920.0 + result * t;
              result = -1216962637349947763.0/17602532858630761920.0 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -89832164536836563.0/82145153340276888960.0;
              result = 3201226073650619.0/488959246073076720.0 + result * t;
              result = -6563969015898899.0/488959246073076720.0 + result * t;
              result = 61800758658685.0/16298641535769224.0 + result * t;
              result = 1981840390203751.0/73343886910961508.0 + result * t;
              result = -854080002959217.0/20373301919711530.0 + result * t;
              result = 1839757479402751.0/91679858638701885.0 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 13853391933823661.0/82145153340276888960.0;
              result = -13002738769221707.0/11735021905753841280.0 + result * t;
              result = 3771644957067319.0/1303891322861537920.0 + result * t;
              result = -7979834368714931.0/2347004381150768256.0 + result * t;
              result = 159017249377589.0/260778264572307584.0 + result * t;
              result = 11301637362262693.0/3911673968584613760.0 + result * t;
              result = -11891175216555961.0/3911673968584613760.0 + result * t;
              result = 81136211446727869.0/82145153340276888960.0 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -122486340027719.0/11735021905753841280.0;
              result = 47258509144553.0/651945661430768960.0 + result * t;
              result = -6751215592079.0/31545757811166240.0 + result * t;
              result = 33756077960395.0/97791849214615344.0 + result * t;
              result = -47258509144553.0/146687773821923016.0 + result * t;
              result = 6751215592079.0/40746603839423060.0 + result * t;
              result = -6751215592079.0/183359717277403770.0 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 964459370297.0/11735021905753841280.0;
              result = -6751215592079.0/11735021905753841280.0 + result * t;
              result = 6751215592079.0/3911673968584613760.0 + result * t;
              result = -6751215592079.0/2347004381150768256.0 + result * t;
              result = 6751215592079.0/2347004381150768256.0 + result * t;
              result = -6751215592079.0/3911673968584613760.0 + result * t;
              result = 6751215592079.0/11735021905753841280.0 + result * t;
              result = -964459370297.0/11735021905753841280.0 + result * t;
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
    return lagrangeSplineBasis.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  LagrangeSplineBasis<LT, IT> lagrangeSplineBasis;
};

// default type-def (unsigned int for level and index)
typedef LagrangeNotAKnotSplineBasis<unsigned int, unsigned int> SLagrangeNotAKnotSplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* LAGRANGE_NOTAKNOT_SPLINE_BASE_HPP */
