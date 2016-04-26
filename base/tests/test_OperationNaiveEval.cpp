// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>
#include <random>

#include "BasisEval.hpp"

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridIndex;
using sgpp::base::GridType;
using sgpp::base::OperationNaiveEval;
using sgpp::base::OperationNaiveEvalGradient;
using sgpp::base::OperationNaiveEvalHessian;
using sgpp::base::OperationNaiveEvalPartialDerivative;
using sgpp::base::SBasis;
using sgpp::base::SPolyBase;
using sgpp::base::SPolyBoundaryBase;

double basisEval(SBasis& basis, GridIndex::level_type l, GridIndex::index_type i,
                        double x) {
  SPolyBase* polyBasis = dynamic_cast<SPolyBase*>(&basis);
  SPolyBoundaryBase* polyBoundaryBasis = dynamic_cast<SPolyBoundaryBase*>(&basis);

  if (polyBasis != nullptr) {
    return polyBasis->evalSave(l, i, x);
  } else if (polyBoundaryBasis != nullptr) {
    return polyBoundaryBasis->evalSave(l, i, x);
  } else {
    return basis.eval(l, i, x);
  }
}

BOOST_AUTO_TEST_CASE(TestOperationNaiveEval) {
  const size_t d = 2;
  const size_t l = 4;
  const size_t p = 3;

  std::mt19937 generator;
  generator.seed(42);
  std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
  std::normal_distribution<double> normalDistribution(0.0, 1.0);

  // Test All The Grids!
  std::vector<std::unique_ptr<Grid>> grids;
  grids.push_back(Grid::createBsplineGrid(d, p));
  grids.push_back(Grid::createBsplineBoundaryGrid(d, p));
  grids.push_back(Grid::createBsplineClenshawCurtisGrid(d, p));
  grids.push_back(Grid::createModBsplineGrid(d, p));
  grids.push_back(Grid::createModBsplineClenshawCurtisGrid(d, p));
  grids.push_back(Grid::createFundamentalSplineGrid(d, p));
  grids.push_back(Grid::createModFundamentalSplineGrid(d, p));
  grids.push_back(Grid::createLinearGrid(d));
  grids.push_back(Grid::createLinearBoundaryGrid(d));
  grids.push_back(Grid::createLinearClenshawCurtisGrid(d));
  grids.push_back(Grid::createModLinearGrid(d));
  grids.push_back(Grid::createWaveletGrid(d));
  grids.push_back(Grid::createWaveletBoundaryGrid(d));
  grids.push_back(Grid::createModWaveletGrid(d));
  grids.push_back(Grid::createPolyGrid(d, p));
  grids.push_back(Grid::createPolyBoundaryGrid(d, p));

  std::vector<std::unique_ptr<SBasis>> bases;
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineBoundaryBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineClenshawCurtisBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SBsplineModifiedClenshawCurtisBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SFundamentalSplineBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SFundamentalSplineModifiedBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearBoundaryBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearClenshawCurtisBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SLinearModifiedBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWaveletBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWaveletBoundaryBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SWaveletModifiedBase()));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyBase(p)));
  bases.push_back(std::unique_ptr<SBasis>(new sgpp::base::SPolyBoundaryBase(p)));

  for (size_t k = 0; k < grids.size(); k++) {
    Grid& grid = *grids[k];
    SBasis& basis = *bases[k];

    // don't test gradients for linear function
    const bool hasGradients = (grid.getType() == GridType::Bspline) ||
                              (grid.getType() == GridType::BsplineBoundary) ||
                              (grid.getType() == GridType::BsplineClenshawCurtis) ||
                              (grid.getType() == GridType::ModBspline) ||
                              (grid.getType() == GridType::ModBsplineClenshawCurtis) ||
                              (grid.getType() == GridType::FundamentalSpline) ||
                              (grid.getType() == GridType::ModFundamentalSpline) ||
                              (grid.getType() == GridType::Wavelet) ||
                              (grid.getType() == GridType::WaveletBoundary) ||
                              (grid.getType() == GridType::ModWavelet);

    // create regular sparse grid
    grid.getGenerator().regular(l);
    const size_t n = grid.getSize();
    DataVector alpha(n);

    for (size_t i = 0; i < n; i++) {
      GridIndex& gp = *grid.getStorage().get(i);

      // don't forget to set the point distribution to Clenshaw-Curtis
      // if necessary (currently not done automatically)
      if ((grid.getType() == GridType::BsplineClenshawCurtis) ||
          (grid.getType() == GridType::ModBsplineClenshawCurtis) ||
          (grid.getType() == GridType::LinearClenshawCurtis)) {
        gp.setPointDistribution(GridIndex::PointDistribution::ClenshawCurtis);
      }

      alpha[i] = normalDistribution(generator);
    }

    // create operations
    std::unique_ptr<OperationNaiveEval> opEval(sgpp::op_factory::createOperationNaiveEval(grid));
    std::unique_ptr<OperationNaiveEvalGradient> opEvalGradient(nullptr);
    std::unique_ptr<OperationNaiveEvalHessian> opEvalHessian(nullptr);
    std::unique_ptr<OperationNaiveEvalPartialDerivative> opEvalPartialDerivative(nullptr);

    if (hasGradients) {
      opEvalGradient = sgpp::op_factory::createOperationNaiveEvalGradient(grid);
      opEvalHessian = sgpp::op_factory::createOperationNaiveEvalHessian(grid);
      opEvalPartialDerivative =
          sgpp::op_factory::createOperationNaiveEvalPartialDerivative(grid);
    }

    DataVector x(d);
    DataVector fxGradient(d);
    DataMatrix fxHessian(d, d);

    for (size_t k = 0; k < 20; k++) {
      // evaluate at random point
      for (size_t t = 0; t < d; t++) {
        x[t] = uniformDistribution(generator);
      }

      double fx = 0.0;
      fxGradient.setAll(0.0);
      fxHessian.setAll(0.0);

      for (size_t i = 0; i < n; i++) {
        // evaluate function by hand
        GridIndex& gp = *grid.getStorage().get(i);
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
              val *= basisEvalDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
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
                val *= basisEvalDxDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              } else if ((t == j) || (t == k)) {
                val *= basisEvalDx(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              } else {
                val *= basisEval(basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              }
            }

            fxHessian(j, k) += val;
          }
        }
      }

      // test function evaluation
      double fx2 = opEval->eval(alpha, x);
      BOOST_CHECK_CLOSE(fx, fx2, 1e-9);

      if (!hasGradients) {
        continue;
      }

      DataVector fxGradient2(d);
      fx2 = opEvalGradient->evalGradient(alpha, x, fxGradient2);

// test function evaluation
      BOOST_CHECK_CLOSE(fx, fx2, 1e-9);

      for (size_t t = 0; t < d; t++) {
        // test gradient evaluation
        BOOST_CHECK_CLOSE(fxGradient[t], fxGradient2[t], 1e-9);

// test partial derivative evaluation
        BOOST_CHECK_CLOSE(opEvalPartialDerivative->evalPartialDerivative(alpha, x, t),
                          fxGradient[t], 1e-9);
      }

      fxGradient2.setAll(0.0);
      DataMatrix fxHessian2(d, d);
      fx2 = opEvalHessian->evalHessian(alpha, x, fxGradient2, fxHessian2);

// test function evaluation
      BOOST_CHECK_CLOSE(fx, fx2, 1e-9);

      for (size_t t1 = 0; t1 < d; t1++) {
        // test gradient evaluation
        BOOST_CHECK_CLOSE(fxGradient[t1], fxGradient2[t1], 1e-9);

        for (size_t t2 = 0; t2 < d; t2++) {
          // test Hessian evaluation
          BOOST_CHECK_CLOSE(fxHessian(t1, t2), fxHessian2(t1, t2), 1e-9);
        }
      }
    }
  }
}
