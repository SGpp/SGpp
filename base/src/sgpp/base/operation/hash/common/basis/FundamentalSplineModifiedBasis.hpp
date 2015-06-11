// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FUNDAMENTAL_SPLINE_MODIFIED_BASE_HPP
#define FUNDAMENTAL_SPLINE_MODIFIED_BASE_HPP

#include <cmath>
#include <stdexcept>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Modified fundamental spline basis on Noboundary grids.
     */
    template <class LT, class IT>
    class FundamentalSplineModifiedBasis: public Basis<LT, IT> {
      public:
        /**
         * Default constructor.
         */
        FundamentalSplineModifiedBasis() :
          fundamentalSplineBasis(FundamentalSplineBasis<LT, IT>()),
          bsplineBasis(BsplineBasis<LT, IT>()),
          bsplineModifiedBasis(BsplineModifiedBasis<LT, IT>()) {
        }

        /**
         * Constructor.
         *
         * @param degree    fundamental spline degree, must be odd
         *                  (if it's even, degree - 1 is used)
         */
        FundamentalSplineModifiedBasis(size_t degree) :
          fundamentalSplineBasis(FundamentalSplineBasis<LT, IT>(degree)),
          bsplineBasis(BsplineBasis<LT, IT>(degree)),
          bsplineModifiedBasis(BsplineModifiedBasis<LT, IT>(degree)) {
          switch (bsplineBasis.getDegree()) {
            case 1:
              coefficients = {2.0, 1.0};
              break;
            case 3:
              coefficients = {
                3.00555349946513493861893415142970724186027337856089,
                1.97927405783630982421797926692905631343500649801563,
                1.07735026918962576450914878085406750439970062937659,
                -0.288675134594812882254574390345326331033809015521984,
                0.0773502691896257645091487805272378197355354327113475,
                -0.0207259421636901757820207317636249479083327153234060,
                0.00555349946513493861893414652726197189779542858227669,
                -0.00148805569684957869371585434542293968284899900570070,
                0.000398723322263376155929270854429786833600567440526131,
                -0.000106837592203925930001229072296207651553270756403819,
                0.0000286270465523275640756454347550437726125155850891456,
                -7.67059400538432630135266672396743889679158395276336e-6,
                2.05532946920974112976523214082598297465075072190782e-6,
                -5.50723871454638217708261839336493001811418934867907e-7,
                1.47566016608811741067815216519989032594925017563811e-7,
                -3.95401949806087465629990267434631285682811353873362e-8,
                1.05947633136232451841808904538634816781995239855340e-8,
                -2.83885827388423417372453507199079814451696055479987e-9,
                7.60669781913691510717249834099710899868318233665478e-10,
                -2.03820853770531869144464264408045454956312379862036e-10,
              };
              break;
            case 5:
              coefficients = {
                3.97761884901317601237505294647935366430548260075455,
                3.05197891743048952938320540794415752252378431262053,
                1.87929705844552798159581670361280439862378201234362,
                1.27994107042222447033988597112101657209527171520227,
                -0.641144390588586790499331240618970433570552138826048,
                0.279941070422224459675698460229065627804305772573088,
                -0.120702941554471989045150037862892980356355531983641,
                0.0519789174304894592205288942891007549201467360258633,
                -0.0223811509868238238226774779024188137457479738263700,
                0.00963678523916236633290992658839119412070080078021938,
                -0.00414936272610310304754576339225598748631415243408659,
                0.00178661332089446143322022453862526514721569525153020,
                -0.000769271651848896958686067511113953123346862349404142,
                0.000331229408555179236582857349959762014249696989546598,
                -0.000142619217560355094053909977433812570893130773236695,
                0.0000614083191042623235506050662723499625441246270561106,
                -0.0000264409083131473972258968002893285943458686486444272,
                0.0000113848032745721823521844968668107540520625936592734,
                -4.90201562161384250942984508867883378040615483109085e-6,
                2.11068707776587432371092034633340130288908179531579e-6,
                -9.08809821128470725014271399616415376068694322895240e-7,
                3.91311104180256378441788970861238025604085853523032e-7,
                -1.68488914506487895598684783698373939882378727632567e-7,
                7.25471728461288350678683430094112219691511858871106e-8,
                -3.12370241293437167236230369277850285984726780190644e-8,
                1.34498925068624274589037267977545049340724248700694e-8,
                -5.79119213459962633826859405412307210223840279586775e-9,
                2.49354456347787225366809386840207060055687387654119e-9,
                -1.07365881592873747751709282016317655415450486554061e-9,
                4.62291017335463009429936619151300495389498474829422e-10,
                -1.99051115250416986152072491494923818169957214087854e-10
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
         * @return      value of modified fundamental spline basis function
         */
        inline float_t eval(LT l, IT i, float_t x) {
          if (l == 1) {
            return 1.0;
          }

          const IT hInv = static_cast<IT>(1) << l;

          if ((i != 1) && (i != hInv - 1)) {
            return fundamentalSplineBasis.eval(l, i, x);
          }

          if (i == hInv - 1) {
            // mirror the situation at x = 0.5
            x = 1.0 - x;
            i = 1;
          }

          const size_t p = bsplineBasis.getDegree();
          const float_t t = x * hInv + static_cast<float_t>(p);

          /*float_t y = 0.0;

          for (size_t k = 0; k < coefficients.size(); k++) {
            y += coefficients[k] * bsplineBasis.uniformBSpline(
                   t - static_cast<float_t>(k), p);
          }*/

          int kMin = std::max(static_cast<int>(t) - static_cast<int>(p + 1) + 1,
                              0);
          int kMax = std::min(static_cast<int>(t),
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
         * @return      value of derivative of modified
         *              fundamental spline basis function
         */
        inline float_t evalDx(LT l, IT i, float_t x) {
          if (l == 1) {
            return 0.0;
          }

          const IT hInv = static_cast<IT>(1) << l;

          if ((i != 1) && (i != hInv - 1)) {
            return fundamentalSplineBasis.evalDx(l, i, x);
          }

          const float_t hInvDbl = static_cast<float_t>(hInv);
          // inner derivative
          float_t dxFactor = hInvDbl;

          if (i == hInv - 1) {
            // mirror the situation at x = 0.5
            x = 1.0 - x;
            i = 1;
            dxFactor *= -1.0;
          }

          const size_t p = bsplineBasis.getDegree();
          const float_t t = x * hInv + static_cast<float_t>(p);

          /*float_t y = 0.0;

          for (size_t k = 0; k < coefficients.size(); k++) {
            y += coefficients[k] * bsplineBasis.uniformBSplineDx(
                   t - static_cast<float_t>(k), p);
          }*/

          int kMin = std::max(static_cast<int>(t) - static_cast<int>(p + 1) + 1,
                              0);
          int kMax = std::min(static_cast<int>(t),
                              static_cast<int>(coefficients.size()) - 1);

          float_t y = 0;

          for (int k = kMin; k <= kMax; k++) {
            y += coefficients[std::abs(k)] * bsplineBasis.uniformBSplineDx(
                   t - static_cast<float_t>(k), p);
          }

          return dxFactor * y;
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of modified
         *              fundamental spline basis function
         */
        inline float_t evalDxDx(LT l, IT i, float_t x) {
          if (l == 1) {
            return 0.0;
          }

          const IT hInv = static_cast<IT>(1) << l;

          if ((i != 1) && (i != hInv - 1)) {
            return fundamentalSplineBasis.evalDxDx(l, i, x);
          }

          const float_t hInvDbl = static_cast<float_t>(hInv);
          // inner derivative
          const float_t dxFactor = hInvDbl * hInvDbl;

          if (i == hInv - 1) {
            // mirror the situation at x = 0.5
            x = 1.0 - x;
            i = 1;
          }

          const size_t p = bsplineBasis.getDegree();
          const float_t t = x * hInv + static_cast<float_t>(p);

          /*float_t y = 0.0;

          for (size_t k = 0; k < coefficients.size(); k++) {
            y += coefficients[k] * bsplineBasis.uniformBSplineDxDx(
                   t - static_cast<float_t>(k), p);
          }*/

          int kMin = std::max(static_cast<int>(t) - static_cast<int>(p + 1) + 1,
                              0);
          int kMax = std::min(static_cast<int>(t),
                              static_cast<int>(coefficients.size()) - 1);

          float_t y = 0;

          for (int k = kMin; k <= kMax; k++) {
            y += coefficients[std::abs(k)] * bsplineBasis.uniformBSplineDxDx(
                   t - static_cast<float_t>(k), p);
          }

          return dxFactor * y;
        }

        /**
         * @return      B-spline degree
         */
        inline size_t getDegree() const {
          return bsplineBasis.getDegree();
        }

      protected:
        /// fundamental spline basis for fundamental spline evaluation
        FundamentalSplineBasis<LT, IT> fundamentalSplineBasis;
        BsplineBasis<LT, IT> bsplineBasis;
        BsplineModifiedBasis<LT, IT> bsplineModifiedBasis;
        std::vector<float_t> coefficients;
    };

    // default type-def (unsigned int for level and index)
    typedef FundamentalSplineModifiedBasis<unsigned int, unsigned int>
    SFundamentalSplineModifiedBase;
  }
}

#endif /* FUNDAMENTAL_SPLINE_MODIFIED_BASE_HPP */
