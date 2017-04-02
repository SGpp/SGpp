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
  BOOST_AUTO_TEST_SUITE_END()
}  // namespace pde
}  // namespace sgpp
