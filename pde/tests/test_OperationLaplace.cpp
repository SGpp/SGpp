// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/globaldef.hpp>
namespace sgpp {
namespace pde {
  /*
    For the BsplineBase we can comapare the OperationLaplace to the Explicit implementation.
    For the other Bases the 1D result is equal to
    = sum_k int (phi'_{i_k}(x_k) * phi'_{j_k}(x_k)) dxk
   */
  BOOST_AUTO_TEST_SUITE(testOperationLaplace)

  BOOST_AUTO_TEST_CASE(testOperationLaplaceBspline) {
    const size_t d = 3;
    const size_t l = 3;
    sgpp::base::Grid* grid(sgpp::base::Grid::createBsplineGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLaplaceExplicit(*grid);
      sgpp::base::OperationMatrix* opImplicit =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector alpha(grid->getSize(), 1.0);
    sgpp::base::DataVector resultImplicit(grid->getSize());
    sgpp::base::DataVector resultExplicit(grid->getSize());

    opExplicit->mult(alpha, resultExplicit);
    opImplicit->mult(alpha, resultImplicit);
    for (size_t i = 0; i < grid->getSize(); i++) {
      BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-15);
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplaceBsplineBoundary1D) {
    const size_t resolution = 10000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createBsplineBoundaryGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SBsplineBoundaryBase basis =
      dynamic_cast<sgpp::base::SBsplineBoundaryBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);
      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 7);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplaceBsplineClenshawCurtis1D) {
    const size_t resolution = 20000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createBsplineClenshawCurtisGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SBsplineClenshawCurtisBase basis =
      dynamic_cast<sgpp::base::SBsplineClenshawCurtisBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 10);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplaceModBsplineClenshawCurtis1D) {
    const size_t resolution = 10000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createModBsplineClenshawCurtisGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SBsplineModifiedClenshawCurtisBase basis =
      dynamic_cast<sgpp::base::SBsplineModifiedClenshawCurtisBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);
      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 7);
      }
    }
  }


  BOOST_AUTO_TEST_CASE(testOperationLaplacePoly1D) {
    const size_t resolution = 10000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createPolyGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SPolyBase basis = dynamic_cast<sgpp::base::SPolyBase&>(grid->getBasis());
    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        BOOST_CHECK_CLOSE(result.get(j), approx, 7);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplacePolyBoundary1D) {
    const size_t resolution = 10000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createPolyBoundaryGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SPolyBoundaryBase basis =
      dynamic_cast<sgpp::base::SPolyBoundaryBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 7);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplaceModPoly1D) {
    const size_t resolution = 10000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createModPolyGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SPolyModifiedBase basis =
      dynamic_cast<sgpp::base::SPolyModifiedBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 7);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplacePolyClenshawCurtis1D) {
    const size_t resolution = 20000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createPolyClenshawCurtisGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SPolyClenshawCurtisBase basis =
      dynamic_cast<sgpp::base::SPolyClenshawCurtisBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 10);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(testOperationLaplacePolyClenshawCurtisBoundary1D) {
    const size_t resolution = 20000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createPolyClenshawCurtisBoundaryGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SPolyClenshawCurtisBoundaryBase basis =
      dynamic_cast<sgpp::base::SPolyClenshawCurtisBoundaryBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 10);
      }
    }
  }


  BOOST_AUTO_TEST_CASE(testOperationLaplaceModPolyClenshawCurtis1D) {
    const size_t resolution = 10000;
    const size_t d = 1;
    const size_t l = 4;
    sgpp::base::Grid* grid(sgpp::base::Grid::createModPolyClenshawCurtisGrid(d, 3));
    grid->getGenerator().regular(l);
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::SPolyModifiedClenshawCurtisBase basis =
      dynamic_cast<sgpp::base::SPolyModifiedClenshawCurtisBase&>(grid->getBasis());

    sgpp::base::OperationMatrix* op =
      sgpp::op_factory::createOperationLaplace(*grid);

    sgpp::base::DataVector result(grid->getSize(), 0);
    for (size_t i = 0; i < grid->getSize(); i++) {
      sgpp::base::DataVector alpha(grid->getSize(), 0.0);
      alpha[i] = 1.0;
      op->mult(alpha, result);

      for (size_t j = 0; j < grid->getSize(); j++) {
        double approx = 0.0;
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(0);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(0);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(0);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(0);
        // trapezoidal rule
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c)/static_cast<double>(resolution);
          approx += basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
        }
        approx += 0.5 * basis.evalDx(lik, iik, 0) * basis.evalDx(ljk, ijk, 0);
        approx += 0.5 * basis.evalDx(lik, iik, 1) * basis.evalDx(ljk, ijk, 1);

        approx /= static_cast<double>(resolution);
        if (result.get(j) > 1e-4)
          BOOST_CHECK_CLOSE(result.get(j), approx, 7);
      }
    }
  }

BOOST_AUTO_TEST_SUITE_END()
}  // namespace pde
}  // namespace sgpp
