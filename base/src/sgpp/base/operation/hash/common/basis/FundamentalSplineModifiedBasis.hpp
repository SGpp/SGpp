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
              coefficients = {1.0};
              break;
            case 3:
              coefficients = {
                1.04674578112205663375656322699728967679562834439521,
                -0.280474686732339802539379361983738060773770066371254,
                0.0751529658073025764009542209376625662994519210898056,
                -0.0201371764968705030644375217669122044240376179879689,
                0.00539574018017943585679586612998625139669855086207003,
                -0.00144578422384724036274594275303280116275658546031119,
                0.000387396715209525594187904882144953254327790979174752,
                -0.000103802636990862014005676775547011854554578456387813,
                0.0000278138327539224618348022200430941638905228463764988,
                -7.45269402482783333353210462536480100751292911818252e-6,
                1.99694334538887149932619845836504013952887009623127e-6,
                -5.35079356727652663772689208095359550602551266742578e-7,
                1.43374081521739155764558374016398062881334970739038e-7,
                -3.84169693593039592855442879702327009227886162135745e-8,
                1.02937959154766813776187778645327408098194941152597e-8,
                -2.75821430260276622493082348789826231648936024746435e-9,
                7.39061294934383522104516087060308456137946874597714e-10,
                -1.98030877134767863487240860342971508062427250926503e-10
              };
              break;
            case 5:
              coefficients = {
                1.11169781137737472588093757238415436343418026149603,
                -0.567837492469660850396483303052716778311637892703584,
                0.248339627548840278715119620887958260566773569273942,
                -0.107094531838723445091464168363929582325059204574711,
                0.0461194024171585616078017058434058157128886956594521,
                -0.0198581852905616219628360799433967715056734936337255,
                0.00855045828010639858177424943774310198162267062157634,
                -0.00368161711309847029032939268175535835778573890614490,
                0.00158521359098153057740011951950482994995508935314007,
                -0.000682553893225497013957595949195807101489497903337790,
                0.000293890879535733221279566679668423262497589117728893,
                -0.000126542167467585554141115849825956467677039716154020,
                0.0000544859376802235887149863969450875235340069835400315,
                -0.0000234603015287334725871379542423464686473794634889565,
                0.0000101014274738046184659300277316765286928486393464088,
                -4.34942564073864181160621177504195949416188862229763e-6,
                1.87275545494656690098008921636071441805797448959485e-6,
                -8.06362330046987386515109253852525621507847175494985e-7,
                3.47199740148325220538430739772932722738175082939161e-7,
                -1.49495648627385837878096330166175147020086680917623e-7,
                6.43691407976723777068271431119458629487503031483672e-8,
                -2.77157651414848679072379101867579030264919222176005e-8,
                1.19337251959362012296680640148996238553485907844255e-8,
                -5.13836786843593259307698280578591013801769497335238e-9,
                2.21245452847915414326784971081810487549275204432462e-9,
                -9.52628376542820646014981500032974866215123088645845e-10,
                4.10178293887209531044927022396464472183043902103230e-10,
                -1.76612661263360302066433255616551731354875783146013e-10
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

          const float_t hInvDbl = static_cast<float_t>(hInv);

          if (i == hInv - 1) {
            // mirror the situation at x = 0.5
            x = 1.0 - x;
            i = 1;
          }

          const size_t p = bsplineBasis.getDegree();
          const float_t t = x * hInv - 1.0 +
                            static_cast<float_t>(p + 1) / 2.0;
          float_t y = coefficients[0] *
                      bsplineModifiedBasis.modifiedBSpline(x * hInvDbl, p);

          for (size_t k = 1; k < coefficients.size(); k++) {
            y += coefficients[k] *
                 bsplineBasis.uniformBSpline(t - static_cast<float_t>(k), p);
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
            return 1.0;
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
          const float_t t = x * hInv - 1.0 +
                            static_cast<float_t>(p + 1) / 2.0;
          float_t y = coefficients[0] *
                      bsplineModifiedBasis.modifiedBSplineDx(x * hInvDbl, p);

          for (size_t k = 1; k < coefficients.size(); k++) {
            y += coefficients[k] *
                 bsplineBasis.uniformBSplineDx(t - static_cast<float_t>(k), p);
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
            return 1.0;
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
          const float_t t = x * hInv - 1.0 +
                            static_cast<float_t>(p + 1) / 2.0;
          float_t y = coefficients[0] *
                      bsplineModifiedBasis.modifiedBSplineDxDx(x * hInvDbl, p);

          for (size_t k = 1; k < coefficients.size(); k++) {
            y += coefficients[k] *
                 bsplineBasis.uniformBSplineDxDx(
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
