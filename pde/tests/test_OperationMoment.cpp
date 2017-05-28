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
    double firstMomentApproximation(sgpp::base::Grid* grid) {
      const size_t d = grid->getDimension();
      const size_t resolution = 10000;
      const double h = 1.0 / static_cast<double>(resolution);
      sgpp::base::GridStorage& storage = grid->getStorage();
      sgpp::base::SBasis& basis = const_cast<sgpp::base::SBasis&>(grid->getBasis());
      double res = 1.0;
      for (size_t k = 0; k < d; k++) {
        const sgpp::base::level_t lik = storage.getPoint(i).getLevel(k);
        const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(k);
        const sgpp::base::index_t iik = storage.getPoint(i).getIndex(k);
        const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(k);
        // trapezoidal rule
        double temp_res = 0.0;
        // --------------------------------------------------------------------------
        // apply trapezoidal rule
        temp_res += x * basis.eval(lik, iik, 0.0) / 2.0;
        for (size_t c = 1; c < resolution; c++) {
          double x = static_cast<double>(c) * h;
          temp_res += x * basis.eval(lik, iik, 0.0);
        }
        temp_res += x * basis.eval(lik, iik, 0.0) / 2.0;
        // --------------------------------------------------------------------------
        res *= temp_res * h;
      }
      return res;
    }

    BOOST_AUTO_TEST_SUITE(testOperationMoment)

    BOOST_AUTO_TEST_CASE(testOperationMomentLinear) {
    }

    BOOST_AUTO_TEST_CASE(testOperationMomentModLinear) {
    }
    BOOST_AUTO_TEST_SUITE_END()
  }  // namespace pde
}  // namespace sgpp
