// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FUNDAMENTAL_SPLINE_BASE_HPP
#define FUNDAMENTAL_SPLINE_BASE_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Fundamental spline basis on Noboundary grids.
     */
    template <class LT, class IT>
    class FundamentalSplineBasis: public Basis<LT, IT> {
      public:
        /**
         * Default constructor.
         */
        FundamentalSplineBasis():
          bsplineBasis(BsplineBasis<LT, IT>()) {
        }

        /**
         * Constructor.
         *
         * @param degree    spline degree, must be odd
         *                  (if it's even, degree - 1 is used)
         */
        FundamentalSplineBasis(size_t degree) :
          bsplineBasis(BsplineBasis<LT, IT>(degree)) {
          switch (bsplineBasis.getDegree()) {
            case 1:
              coefficients = {1.0};
              break;
            case 3:
              coefficients = {
                1.73205080756887729352744634150587236694280525381038,
                -0.464101615137754587054892683011744733885610507620761,
                0.124355652982141054692124390541106568599636776672664,
                -0.0333209967908096317136048791526815405129365990698963,
                0.00892833418109747216229512606961959345210961960692092,
                -0.00239233993358025693557562512579683329550187935778736,
                0.000641025553223555580007374433567739729897897824228503,
                -0.000171762279313965384453872608474125624089711939126657,
                0.0000460235640323059578081160003287627664609499322781231,
                -0.0000123319768152584467785913928409254417540877899858359,
                3.30434322872782930624957103493900055540122766522066e-6,
                -8.85396099652870446406891298830560467517120675046676e-7,
                2.37241169883652479377994160383241314667255034966049e-7,
                -6.35685798817394711050853427024047911518994648175213e-8,
                1.70331496433054050423472104263778499403428243040360e-8,
                -4.56401869148214906430349900310660860947183239862251e-9,
                1.22292512262319121486678558604858449754450529045408e-9,
                -3.27681799010615795163643341087729380706188763193797e-10
              };
              break;
            case 5:
              coefficients = {
                2.84217092202162248422242525148418205475658661917667,
                -1.32172947298750768735998028776961192189206388349626,
                0.573325870961657892019454183031902162226302538072862,
                -0.247041927402274729170186016854190874786837894238021,
                0.106378004643299472277922588629147951958892089105880,
                -0.0458040841912516591395386630410933109916954034416355,
                0.0197221240122630336904966447217218387286853206101436,
                -0.00849186101974092279916228955626515763827813123063995,
                0.00365638603314743455552362726853087632718865210380614,
                -0.00157434968651961051227729799197687633467122589882045,
                0.000677876162780151742373231809533206973546925904621890,
                -0.000291876764082027135854379002927170092545987666338926,
                0.000125674939005129299266675169782670526361790885292975,
                -0.0000541125304839056039712692770427125549202246605969030,
                0.0000232995215955657414930622640274871206963961663733109,
                -0.0000100321995982740617950219913179001062320480775889102,
                4.31961782420307211478838177910869763515332689964936e-6,
                -1.85992094399547136139608694737171038488904363199298e-6,
                8.00836106039360997774573448406067423849409993312779e-7,
                -3.44820284328089279473106070950899899346201916568684e-7,
                1.48471113611678694688789102478719780901997804310541e-7,
                -6.39280012776681974263386928733635817982326126802284e-8,
                2.75258213395395521509096474292775463740326826632236e-8,
                -1.18519400774841082207076970150342111472383078541184e-8,
                5.10315321267081017438030428730477499408907380972312e-9,
                -2.19729196585008046449204876551883805520473724456945e-9,
                9.46099750875882155183659902986583279607618460934836e-10,
                -4.07367228624581691480153408297792325363033041404921e-10,
                1.75402285862183638775398634813118446570135203911292e-10
              };
              break;
            default:
              throw std::invalid_argument("Degree is unsupported.");
          }
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of fundamental spline basis function
         */
        inline float_t eval(LT l, IT i, float_t x) {
          const size_t p = bsplineBasis.getDegree();
          const float_t hInv = static_cast<float_t>(static_cast<IT>(1) << l);
          const float_t t = x * hInv - static_cast<float_t>(i) +
                            static_cast<float_t>(p + 1) / 2.0;

          /*float_t y = coefficients[0] * bsplineBasis.uniformBSpline(t, p);

          for (size_t k = 1; k < coefficients.size(); k++) {
            y += coefficients[k] * (
                   bsplineBasis.uniformBSpline(
                     t - static_cast<float_t>(k), p) +
                   bsplineBasis.uniformBSpline(
                     t + static_cast<float_t>(k), p));
          }*/

          int kMin = std::max(static_cast<int>(std::floor(t)) -
                              static_cast<int>(p + 1) + 1,
                              1 - static_cast<int>(coefficients.size()));
          int kMax = std::min(static_cast<int>(std::floor(t)),
                              static_cast<int>(coefficients.size()) - 1);

          float_t y = 0;

          for (int k = kMin; k <= kMax; k++) {
            y += coefficients[std::abs(k)] * bsplineBasis.uniformBSpline(
                   t - static_cast<float_t>(k), p);
          }

          return y;
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of derivative of fundamental spline basis
         *              function
         */
        inline float_t evalDx(LT l, IT i, float_t x) {
          const size_t p = bsplineBasis.getDegree();
          const float_t hInv = static_cast<float_t>(static_cast<IT>(1) << l);
          const float_t t = x * hInv - static_cast<float_t>(i) +
                            static_cast<float_t>(p + 1) / 2.0;

          /*float_t y = coefficients[0] *
                      bsplineBasis.uniformBSplineDx(t, p);

          for (size_t k = 1; k < coefficients.size(); k++) {
            y += coefficients[k] * (
                   bsplineBasis.uniformBSplineDx(
                     t - static_cast<float_t>(k), p) +
                   bsplineBasis.uniformBSplineDx(
                     t + static_cast<float_t>(k), p));
          }*/

          int kMin = std::max(static_cast<int>(std::floor(t)) -
                              static_cast<int>(p + 1) + 1,
                              1 - static_cast<int>(coefficients.size()));
          int kMax = std::min(static_cast<int>(std::floor(t)),
                              static_cast<int>(coefficients.size()) - 1);

          float_t y = 0;

          for (int k = kMin; k <= kMax; k++) {
            y += coefficients[std::abs(k)] * bsplineBasis.uniformBSplineDx(
                   t - static_cast<float_t>(k), p);
          }

          // don't forget the inner derivative
          return hInv * y;
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of fundamental spline basis
         *              function
         */
        inline float_t evalDxDx(LT l, IT i, float_t x) {
          const size_t p = bsplineBasis.getDegree();
          const float_t hInv = static_cast<float_t>(static_cast<IT>(1) << l);
          const float_t t = x * hInv - static_cast<float_t>(i) +
                            static_cast<float_t>(p + 1) / 2.0;

          /*float_t y = coefficients[0] *
                      bsplineBasis.uniformBSplineDxDx(t, p);

          for (size_t k = 1; k < coefficients.size(); k++) {
            y += coefficients[k] * (
                   bsplineBasis.uniformBSplineDxDx(
                     t - static_cast<float_t>(k), p) +
                   bsplineBasis.uniformBSplineDxDx(
                     t + static_cast<float_t>(k), p));
          }*/

          int kMin = std::max(static_cast<int>(std::floor(t)) -
                              static_cast<int>(p + 1) + 1,
                              1 - static_cast<int>(coefficients.size()));
          int kMax = std::min(static_cast<int>(std::floor(t)),
                              static_cast<int>(coefficients.size()) - 1);

          float_t y = 0;

          for (int k = kMin; k <= kMax; k++) {
            y += coefficients[std::abs(k)] * bsplineBasis.uniformBSplineDxDx(
                   t - static_cast<float_t>(k), p);
          }

          // don't forget the inner derivative
          return hInv * hInv * y;
        }

        /**
         * @return      fundamental spline degree
         */
        inline size_t getDegree() const {
          return bsplineBasis.getDegree();
        }

      protected:
        std::vector<float_t> coefficients;
        BsplineBasis<LT, IT> bsplineBasis;
    };

    // default type-def (unsigned int for level and index)
    typedef FundamentalSplineBasis<unsigned int, unsigned int>
    SFundamentalSplineBase;
  }
}

#endif /* FUNDAMENTAL_SPLINE_BASE_HPP */
