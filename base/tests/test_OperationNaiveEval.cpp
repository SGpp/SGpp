#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>
#include <random>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

using namespace SGPP;
using namespace SGPP::base;

SGPP::float_t basisEval(SBasis& basis,
                        GridIndex::level_type l,
                        GridIndex::index_type i,
                        SGPP::float_t x) {
  SPolyBase* polyBasis =
    dynamic_cast<SPolyBase*>(&basis);
  SPolyBoundaryBase* polyBoundaryBasis =
    dynamic_cast<SPolyBoundaryBase*>(&basis);

  if (polyBasis != nullptr) {
    return polyBasis->evalSave(l, i, x);
  } else if (polyBoundaryBasis != nullptr) {
    return polyBoundaryBasis->evalSave(l, i, x);
  } else {
    return basis.eval(l, i, x);
  }
}

SGPP::float_t basisEvalDx(SBasis& basis,
                          GridIndex::level_type l,
                          GridIndex::index_type i,
                          SGPP::float_t x) {
  SBsplineBase* bsplineBasis =
    dynamic_cast<SBsplineBase*>(&basis);
  SBsplineBoundaryBase* bsplineBoundaryBasis =
    dynamic_cast<SBsplineBoundaryBase*>(&basis);
  SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
    dynamic_cast<SBsplineClenshawCurtisBase*>(&basis);
  SBsplineModifiedBase* bsplineModifiedBasis =
    dynamic_cast<SBsplineModifiedBase*>(&basis);
  SFundamentalSplineBase* fundamentalSplineBasis =
    dynamic_cast<SFundamentalSplineBase*>(&basis);
  SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
    dynamic_cast<SFundamentalSplineModifiedBase*>(&basis);
  SWaveletBase* waveletBasis =
    dynamic_cast<SWaveletBase*>(&basis);
  SWaveletBoundaryBase* waveletBoundaryBasis =
    dynamic_cast<SWaveletBoundaryBase*>(&basis);
  SWaveletModifiedBase* waveletModifiedBasis =
    dynamic_cast<SWaveletModifiedBase*>(&basis);

  if (bsplineBasis != nullptr) {
    return bsplineBasis->evalDx(l, i, x);
  } else if (bsplineBoundaryBasis != nullptr) {
    return bsplineBoundaryBasis->evalDx(l, i, x);
  } else if (bsplineClenshawCurtisBasis != nullptr) {
    return bsplineClenshawCurtisBasis->evalDx(l, i, x);
  } else if (bsplineModifiedBasis != nullptr) {
    return bsplineModifiedBasis->evalDx(l, i, x);
  } else if (fundamentalSplineBasis != nullptr) {
    return fundamentalSplineBasis->evalDx(l, i, x);
  } else if (fundamentalSplineModifiedBasis != nullptr) {
    return fundamentalSplineModifiedBasis->evalDx(l, i, x);
  } else if (waveletBasis != nullptr) {
    return waveletBasis->evalDx(l, i, x);
  } else if (waveletBoundaryBasis != nullptr) {
    return waveletBoundaryBasis->evalDx(l, i, x);
  } else if (waveletModifiedBasis != nullptr) {
    return waveletModifiedBasis->evalDx(l, i, x);
  }

  return NAN;
}

SGPP::float_t basisEvalDxDx(SBasis& basis,
                            GridIndex::level_type l,
                            GridIndex::index_type i,
                            SGPP::float_t x) {
  SBsplineBase* bsplineBasis =
    dynamic_cast<SBsplineBase*>(&basis);
  SBsplineBoundaryBase* bsplineBoundaryBasis =
    dynamic_cast<SBsplineBoundaryBase*>(&basis);
  SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
    dynamic_cast<SBsplineClenshawCurtisBase*>(&basis);
  SBsplineModifiedBase* bsplineModifiedBasis =
    dynamic_cast<SBsplineModifiedBase*>(&basis);
  SFundamentalSplineBase* fundamentalSplineBasis =
    dynamic_cast<SFundamentalSplineBase*>(&basis);
  SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
    dynamic_cast<SFundamentalSplineModifiedBase*>(&basis);
  SWaveletBase* waveletBasis =
    dynamic_cast<SWaveletBase*>(&basis);
  SWaveletBoundaryBase* waveletBoundaryBasis =
    dynamic_cast<SWaveletBoundaryBase*>(&basis);
  SWaveletModifiedBase* waveletModifiedBasis =
    dynamic_cast<SWaveletModifiedBase*>(&basis);

  if (bsplineBasis != nullptr) {
    return bsplineBasis->evalDxDx(l, i, x);
  } else if (bsplineBoundaryBasis != nullptr) {
    return bsplineBoundaryBasis->evalDxDx(l, i, x);
  } else if (bsplineClenshawCurtisBasis != nullptr) {
    return bsplineClenshawCurtisBasis->evalDxDx(l, i, x);
  } else if (bsplineModifiedBasis != nullptr) {
    return bsplineModifiedBasis->evalDxDx(l, i, x);
  } else if (fundamentalSplineBasis != nullptr) {
    return fundamentalSplineBasis->evalDxDx(l, i, x);
  } else if (fundamentalSplineModifiedBasis != nullptr) {
    return fundamentalSplineModifiedBasis->evalDxDx(l, i, x);
  } else if (waveletBasis != nullptr) {
    return waveletBasis->evalDxDx(l, i, x);
  } else if (waveletBoundaryBasis != nullptr) {
    return waveletBoundaryBasis->evalDxDx(l, i, x);
  } else if (waveletModifiedBasis != nullptr) {
    return waveletModifiedBasis->evalDxDx(l, i, x);
  }

  return NAN;
}

BOOST_AUTO_TEST_CASE(TestOperationNaiveEval) {
  const size_t d = 2;
  const size_t l = 4;
  const size_t p = 3;

  std::mt19937 generator;
  generator.seed(42);
  std::uniform_real_distribution<SGPP::float_t> uniformDistribution(0.0, 1.0);
  std::normal_distribution<SGPP::float_t> normalDistribution(0.0, 1.0);

  // Test All The Grids!
  std::vector<std::unique_ptr<Grid>> grids;
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createBsplineTruncatedBoundaryGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createModBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createModFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createLinearTruncatedBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createLinearClenshawCurtisGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createModLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createWaveletGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createWaveletTruncatedBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createModWaveletGrid(d))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createPolyGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<Grid>(
                              Grid::createPolyTruncatedBoundaryGrid(d, p))));

  std::vector<std::unique_ptr<SBasis>> bases;
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SBsplineBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SBsplineBoundaryBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SBsplineClenshawCurtisBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SBsplineModifiedBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SFundamentalSplineBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SFundamentalSplineModifiedBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SLinearBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SLinearBoundaryBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SLinearClenshawCurtisBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SLinearModifiedBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SWaveletBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SWaveletBoundaryBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SWaveletModifiedBase())));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SPolyBase(p))));
  bases.push_back(std::move(std::unique_ptr<SBasis>(
                              new SPolyBoundaryBase(p))));

  for (size_t k = 0; k < grids.size(); k++) {
    Grid& grid = *grids[k];
    SBasis& basis = *bases[k];

    // don't test gradients for linear function
    const bool hasGradients =
      (std::string(grid.getType()).find("linear") == std::string::npos) &&
      (std::string(grid.getType()).find("poly") == std::string::npos);

    // create regular sparse grid
    std::unique_ptr<GridGenerator> gridGen(grid.createGridGenerator());
    gridGen->regular(l);
    const size_t n = grid.getSize();
    DataVector alpha(n);

    for (size_t i = 0; i < n; i++) {
      GridIndex& gp = *grid.getStorage()->get(i);

      // don't forget to set the point distribution to Clenshaw-Curtis
      // if necessary (currently not done automatically)
      if ((std::string(grid.getType()) == "bsplineClenshawCurtis") ||
          (std::string(grid.getType()) == "linearClenshawCurtis")) {
        gp.setPointDistribution(GridIndex::PointDistribution::ClenshawCurtis);
      }

      alpha[i] = normalDistribution(generator);
    }

    // create operations
    std::unique_ptr<OperationNaiveEval> opEval(
      op_factory::createOperationNaiveEval(grid));
    std::unique_ptr<OperationNaiveEvalGradient> opEvalGradient(nullptr);
    std::unique_ptr<OperationNaiveEvalHessian> opEvalHessian(nullptr);
    std::unique_ptr<OperationNaiveEvalPartialDerivative>
    opEvalPartialDerivative(nullptr);

    if (hasGradients) {
      opEvalGradient.reset(op_factory::createOperationNaiveEvalGradient(grid));
      opEvalHessian.reset(op_factory::createOperationNaiveEvalHessian(grid));
      opEvalPartialDerivative.reset(
        op_factory::createOperationNaiveEvalPartialDerivative(grid));
    }

    base::DataVector x(d);
    base::DataVector fxGradient(d);
    base::DataMatrix fxHessian(d, d);

    for (size_t k = 0; k < 20; k++) {
      // evaluate at random point
      for (size_t t = 0; t < d; t++) {
        x[t] = uniformDistribution(generator);
      }

      SGPP::float_t fx = 0.0;
      fxGradient.setAll(0.0);
      fxHessian.setAll(0.0);

      for (size_t i = 0; i < n; i++) {
        // evaluate function by hand
        GridIndex& gp = *grid.getStorage()->get(i);
        SGPP::float_t val = alpha[i];

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
                val *= basisEvalDxDx(
                         basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              } else if ((t == j) || (t == k)) {
                val *= basisEvalDx(
                         basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              } else {
                val *= basisEval(
                         basis, gp.getLevel(t), gp.getIndex(t), x[t]);
              }
            }

            fxHessian.set(j, k, fxHessian.get(j, k) + val);
          }
        }
      }

      // test function evaluation
      SGPP::float_t fx2 = opEval->eval(alpha, x);
      BOOST_CHECK_CLOSE(fx, fx2, 1e-10);

      if (!hasGradients) {
        continue;
      }

      base::DataVector fxGradient2(d);
      fx2 = opEvalGradient->evalGradient(alpha, x, fxGradient2);

      // test function evaluation
      BOOST_CHECK_CLOSE(fx, fx2, 1e-10);

      for (size_t t = 0; t < d; t++) {
        // test gradient evaluation
        BOOST_CHECK_CLOSE(fxGradient[t], fxGradient2[t], 1e-10);

        // test partial derivative evaluation
        BOOST_CHECK_CLOSE(opEvalPartialDerivative->evalPartialDerivative(
                            alpha, x, t), fxGradient[t], 1e-10);
      }

      fxGradient2.setAll(0.0);
      base::DataMatrix fxHessian2(d, d);
      fx2 = opEvalHessian->evalHessian(alpha, x, fxGradient2, fxHessian2);

      // test function evaluation
      BOOST_CHECK_CLOSE(fx, fx2, 1e-10);

      for (size_t t1 = 0; t1 < d; t1++) {
        // test gradient evaluation
        BOOST_CHECK_CLOSE(fxGradient[t1], fxGradient2[t1], 1e-10);

        for (size_t t2 = 0; t2 < d; t2++) {
          // test Hessian evaluation
          BOOST_CHECK_CLOSE(fxHessian.get(t1, t2),
                            fxHessian2.get(t1, t2), 1e-10);
        }
      }
    }
  }
}
