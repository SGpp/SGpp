// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBasis.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>
#include <random>

#include "BasisEval.hpp"

using sgpp::base::BoundingBox;
using sgpp::base::BoundingBox1D;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridType;
using sgpp::base::OperationEval;
using sgpp::base::OperationEvalGradient;
using sgpp::base::OperationEvalHessian;
using sgpp::base::OperationEvalPartialDerivative;
using sgpp::base::SBasis;
using sgpp::base::SPolyBase;
using sgpp::base::SPolyBoundaryBase;
using sgpp::base::SPolyModifiedBase;
using sgpp::base::SPolyClenshawCurtisBoundaryBase;

double basisEval(SBasis& basis, GridPoint::level_type l, GridPoint::index_type i, double x) {
  return basis.eval(l, i, x);
}

void checkClose(double x, double y, double tol = 1e-8) { BOOST_CHECK_CLOSE(x, y, tol); }

void checkClose(const DataVector& x, const DataVector& y, double tol = 1e-8) {
  BOOST_CHECK_EQUAL(x.getSize(), y.getSize());

  for (size_t i = 0; i < x.getSize(); i++) {
    BOOST_CHECK_CLOSE(x[i], y[i], tol);
  }
}

void checkClose(const DataMatrix& x, const DataMatrix& y, double tol = 1e-8) {
  BOOST_CHECK_EQUAL(x.getNrows(), y.getNrows());
  BOOST_CHECK_EQUAL(x.getNcols(), y.getNcols());

  for (size_t i = 0; i < x.getNrows(); i++) {
    for (size_t j = 0; j < x.getNcols(); j++) {
      BOOST_CHECK_CLOSE(x(i, j), y(i, j), tol);
    }
  }
}

void checkClose(const std::vector<DataMatrix>& x, const std::vector<DataMatrix>& y,
                double tol = 5e-8) {
  BOOST_CHECK_EQUAL(x.size(), y.size());

  for (size_t k = 0; k < x.size(); k++) {
    BOOST_CHECK_EQUAL(x[k].getNrows(), y[k].getNrows());
    BOOST_CHECK_EQUAL(x[k].getNcols(), y[k].getNcols());

    for (size_t i = 0; i < x[k].getNrows(); i++) {
      for (size_t j = 0; j < x[k].getNcols(); j++) {
        BOOST_CHECK_CLOSE(x[k](i, j), y[k](i, j), tol);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestOperationEvalNaive) {
  const size_t d = 2;
  const size_t l = 4;
  const size_t p = 3;
  const size_t m = 3;
  const size_t N = 20;

  std::mt19937 generator;
  generator.seed(42);
  std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
  std::normal_distribution<double> normalDistribution(0.0, 1.0);

  // Test All The Grids!
  std::vector<std::unique_ptr<Grid>> grids;
  grids.push_back(std::unique_ptr<Grid>(Grid::createBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createBsplineClenshawCurtisGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModBsplineClenshawCurtisGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createFundamentalNakSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createFundamentalSplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createFundamentalSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModFundamentalSplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createWeaklyFundamentalNakSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModWeaklyFundamentalNakSplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createWeaklyFundamentalSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createLinearGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createLinearBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createLinearClenshawCurtisGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createLinearClenshawCurtisBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModLinearGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createNaturalBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createNakBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModNakBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createWaveletGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createWaveletBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModWaveletGrid(d)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createPolyGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createPolyBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModPolyGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createPolyClenshawCurtisBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createPolyClenshawCurtisGrid(d, p)));

  std::vector<std::unique_ptr<SBasis>> bases;
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineBoundaryBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineClenshawCurtisBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineModifiedClenshawCurtisBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SFundamentalNakSplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SFundamentalSplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SFundamentalSplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SFundamentalSplineModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWeaklyFundamentalNakSplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWeaklyFundamentalNakSplineModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWeaklyFundamentalSplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearBoundaryBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearClenshawCurtisBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearClenshawCurtisBoundaryBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearModifiedBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SNaturalBsplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SNakBsplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SNakBsplineModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWaveletBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWaveletBoundaryBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWaveletModifiedBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyBoundaryBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyClenshawCurtisBoundaryBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyClenshawCurtisBase(p)));

  for (size_t k = 0; k < grids.size(); k++) {
    Grid& grid = *grids[k];
    SBasis& basis = *bases[k];

    // only test gradients for bases with derivatives
    const bool hasGradients =
        (grid.getType() == GridType::Bspline) ||
        (grid.getType() == GridType::BsplineBoundary) ||
        (grid.getType() == GridType::BsplineClenshawCurtis) ||
        (grid.getType() == GridType::ModBspline) ||
        (grid.getType() == GridType::ModBsplineClenshawCurtis) ||
        (grid.getType() == GridType::FundamentalNakSplineBoundary) ||
        (grid.getType() == GridType::FundamentalSpline) ||
        (grid.getType() == GridType::FundamentalSplineBoundary) ||
        (grid.getType() == GridType::ModFundamentalSpline) ||
        (grid.getType() == GridType::WeaklyFundamentalNakSplineBoundary) ||
        (grid.getType() == GridType::ModWeaklyFundamentalNakSpline) ||
        (grid.getType() == GridType::WeaklyFundamentalSplineBoundary) ||
        (grid.getType() == GridType::NakBsplineBoundary) ||
        (grid.getType() == GridType::ModNakBspline) ||
        (grid.getType() == GridType::Wavelet) ||
        (grid.getType() == GridType::WaveletBoundary) ||
        (grid.getType() == GridType::ModWavelet);

    // create regular sparse grid
    grid.getGenerator().regular(l);
    const size_t n = grid.getSize();

    // set random bounding box
    BoundingBox& boundingBox = grid.getBoundingBox();
    DataVector innerDerivative(d);

    for (size_t t = 0; t < d; t++) {
      // ensure left < right
      const double left = normalDistribution(generator);
      const double right = left + std::abs(normalDistribution(generator));
      BoundingBox1D boundingBox1D(left, right);
      boundingBox.setBoundary(t, boundingBox1D);
      innerDerivative[t] = 1.0 / (right - left);
    }

    // create operations
    std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEvalNaive(grid));
    std::unique_ptr<OperationEvalGradient> opEvalGradient(nullptr);
    std::unique_ptr<OperationEvalHessian> opEvalHessian(nullptr);
    std::unique_ptr<OperationEvalPartialDerivative> opEvalPartialDerivative(nullptr);

    if (hasGradients) {
      opEvalGradient.reset(sgpp::op_factory::createOperationEvalGradientNaive(grid));
      opEvalHessian.reset(sgpp::op_factory::createOperationEvalHessianNaive(grid));
      opEvalPartialDerivative.reset(
          sgpp::op_factory::createOperationEvalPartialDerivativeNaive(grid));
    }

    // test vector version (single coefficient vector)
    {
      // create coefficient vector
      DataVector alpha(n);

      for (size_t i = 0; i < n; i++) {
        alpha[i] = normalDistribution(generator);
      }

      // x is in unit cube, y is in BoundingBox
      DataVector x(d), y(d);
      double fx;
      DataVector fxGradient(d);
      DataMatrix fxHessian(d, d);

      for (size_t r = 0; r < N; r++) {
        // evaluate at random point
        for (size_t t = 0; t < d; t++) {
          x[t] = uniformDistribution(generator);
          y[t] = boundingBox.getIntervalOffset(t) + boundingBox.getIntervalWidth(t) * x[t];
        }

        fx = 0.0;
        fxGradient.setAll(0.0);
        fxHessian.setAll(0.0);

        for (size_t i = 0; i < n; i++) {
          // evaluate function by hand
          GridPoint& gp = grid.getStorage().getPoint(i);
          double val = alpha[i];

          for (size_t t = 0; t < d; t++) {
            val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
          }

          fx += val;

          if (!hasGradients) {
            continue;
          }

          // evaluate gradient by hand
          for (size_t j = 0; j < d; j++) {
            val = alpha[i];

            for (size_t t = 0; t < d; t++) {
              if (t == j) {
                val *=
                    basisEvalDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]) * innerDerivative[t];
              } else {
                val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              }
            }

            fxGradient[j] += val;
          }

          // evaluate Hessian by hand
          for (size_t j = 0; j < d; j++) {
            for (size_t k = 0; k < d; k++) {
              val = alpha[i];

              for (size_t t = 0; t < d; t++) {
                if ((t == j) && (t == k)) {
                  val *= basisEvalDxDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]) *
                         innerDerivative[t] * innerDerivative[t];
                } else if ((t == j) || (t == k)) {
                  val *=
                      basisEvalDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]) * innerDerivative[t];
                } else {
                  val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
                }
              }

              fxHessian(j, k) += val;
            }
          }
        }

        // test function evaluation
        double fx2 = opEval->eval(alpha, y);
        checkClose(fx, fx2);

        if (!hasGradients) {
          continue;
        }

        // test gradient evaluation
        DataVector fxGradient2(d);
        fx2 = opEvalGradient->evalGradient(alpha, y, fxGradient2);
        checkClose(fx, fx2);
        checkClose(fxGradient, fxGradient2);

        // test partial derivative evaluation
        for (size_t t = 0; t < d; t++) {
          double partDeriv2;
          fx2 = opEvalPartialDerivative->evalPartialDerivative(alpha, y, t, partDeriv2);
          checkClose(fx, fx2);
          checkClose(fxGradient[t], partDeriv2);
        }

        // test Hessian evaluation
        fxGradient2.setAll(0.0);
        DataMatrix fxHessian2(d, d);
        fx2 = opEvalHessian->evalHessian(alpha, y, fxGradient2, fxHessian2);
        checkClose(fx, fx2);
        checkClose(fxGradient, fxGradient2);
        checkClose(fxHessian, fxHessian2);
      }
    }

    // test matrix version (multiple coefficient vectors)
    {
      // create coefficient vector
      DataMatrix alpha(n, m);

      for (size_t i = 0; i < n; i++) {
        for (size_t q = 0; q < m; q++) {
          alpha(i, q) = normalDistribution(generator);
        }
      }

      // x is in unit cube, y is in BoundingBox
      DataVector x(d), y(d);
      DataVector fx(m);
      DataMatrix fxGradient(m, d);
      std::vector<DataMatrix> fxHessian(m, DataMatrix(d, d));

      for (size_t r = 0; r < N; r++) {
        // evaluate at random point
        for (size_t t = 0; t < d; t++) {
          x[t] = uniformDistribution(generator);
          y[t] = boundingBox.getIntervalOffset(t) + boundingBox.getIntervalWidth(t) * x[t];
        }

        fx.setAll(0.0);
        fxGradient.setAll(0.0);

        for (size_t q = 0; q < m; q++) {
          fxHessian[q].setAll(0.0);

          for (size_t i = 0; i < n; i++) {
            // evaluate function by hand
            GridPoint& gp = grid.getStorage().getPoint(i);
            double val = alpha(i, q);

            for (size_t t = 0; t < d; t++) {
              val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
            }

            fx[q] += val;

            if (!hasGradients) {
              continue;
            }

            // evaluate gradient by hand
            for (size_t j = 0; j < d; j++) {
              val = alpha(i, q);

              for (size_t t = 0; t < d; t++) {
                if (t == j) {
                  val *=
                      basisEvalDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]) * innerDerivative[t];
                } else {
                  val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
                }
              }

              fxGradient(q, j) += val;
            }

            // evaluate Hessian by hand
            for (size_t j = 0; j < d; j++) {
              for (size_t k = 0; k < d; k++) {
                val = alpha(i, q);

                for (size_t t = 0; t < d; t++) {
                  if ((t == j) && (t == k)) {
                    val *= basisEvalDxDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]) *
                           innerDerivative[t] * innerDerivative[t];
                  } else if ((t == j) || (t == k)) {
                    val *= basisEvalDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]) *
                           innerDerivative[t];
                  } else {
                    val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
                  }
                }

                fxHessian[q](j, k) += val;
              }
            }
          }
        }

        // test function evaluation
        DataVector fx2(m);
        opEval->eval(alpha, y, fx2);
        checkClose(fx, fx2);

        if (!hasGradients) {
          continue;
        }

        fx2.setAll(0.0);
        DataMatrix fxGradient2(m, d);
        opEvalGradient->evalGradient(alpha, y, fx2, fxGradient2);

        // test gradient evaluation
        checkClose(fx, fx2);
        checkClose(fxGradient, fxGradient2);

        // test partial derivative evaluation
        for (size_t t = 0; t < d; t++) {
          DataVector fxPartDeriv(m);
          fxGradient.getColumn(t, fxPartDeriv);

          fx2.setAll(0.0);
          DataVector fxPartDeriv2(m);
          opEvalPartialDerivative->evalPartialDerivative(alpha, y, t, fx2, fxPartDeriv2);

          checkClose(fx, fx2);
          checkClose(fxPartDeriv, fxPartDeriv2);
        }

        // test Hessian evaluation
        fx2.setAll(0.0);
        fxGradient2.setAll(0.0);
        std::vector<DataMatrix> fxHessian2(m, DataMatrix(d, d));
        opEvalHessian->evalHessian(alpha, y, fx2, fxGradient2, fxHessian2);

        // test function evaluation
        checkClose(fx, fx2);
        checkClose(fxGradient, fxGradient2);
        checkClose(fxHessian, fxHessian2);
      }
    }
  }
}
