// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppCombigridModule

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/combigrid/CombinationGrid.hpp>
#include <sgpp/combigrid/FullGrid.hpp>
#include <sgpp/combigrid/HeterogeneousBasis.hpp>
#include <sgpp/combigrid/IndexVectorRange.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/OperationEvalCombinationGrid.hpp>
#include <sgpp/combigrid/OperationPole.hpp>
#include <sgpp/combigrid/OperationPoleHierarchisationGeneral.hpp>
#include <sgpp/combigrid/OperationPoleHierarchisationLinear.hpp>
#include <sgpp/combigrid/OperationPoleNodalisationBspline.hpp>
#include <sgpp/combigrid/OperationUPFullGrid.hpp>
#include <sgpp/combigrid/OperationUPCombinationGrid.hpp>

#include <boost/test/unit_test.hpp>

#include <memory>
#include <vector>

// TODO(valentjn)
#include <iostream>

using sgpp::base::DataVector;
using sgpp::combigrid::CombinationGrid;
using sgpp::combigrid::FullGrid;
using sgpp::combigrid::HeterogeneousBasis;
using sgpp::combigrid::IndexVector;
using sgpp::combigrid::IndexVectorRange;
using sgpp::combigrid::LevelVector;
using sgpp::combigrid::OperationEvalCombinationGrid;
using sgpp::combigrid::OperationPole;
using sgpp::combigrid::OperationPoleHierarchisationGeneral;
using sgpp::combigrid::OperationPoleHierarchisationLinear;
using sgpp::combigrid::OperationPoleNodalisationBspline;
using sgpp::combigrid::OperationUPFullGrid;
using sgpp::combigrid::OperationUPCombinationGrid;

// fix for clang (from https://stackoverflow.com/a/33755176)
#ifdef __clang__
#include <string>

namespace boost {
namespace unit_test {
namespace ut_detail {

std::string normalize_test_case_name(const_string name) {
  return ((name[0] == '&') ? std::string(name.begin() + 1, name.size() - 1) :
                             std::string(name.begin(), name.size()));
}

}  // namespace ut_detail
}  // namespace unit_test
}  // namespace boost
#endif

BOOST_AUTO_TEST_CASE(testFullGrid) {
  sgpp::base::SBsplineBase basis1d;
  HeterogeneousBasis basis(3, basis1d);
  FullGrid fullGrid({2, 3, 2}, basis);
  IndexVectorRange range(fullGrid);
  BOOST_CHECK_EQUAL(std::vector<IndexVector>(range.begin(), range.end()).size(), 225);

  /*for (const IndexVector& index : IndexVectorRange(fullGrid)) {
    std::cout << index.size() << ": " << index[0] << ", " << index[1] << ", " << index[2] << "\n";
  }*/
}

BOOST_AUTO_TEST_CASE(testCombinationGrid) {
  BOOST_CHECK_EQUAL(CombinationGrid::enumerateLevelsWithSumWithBoundary(3, 5).size(), 21);
  BOOST_CHECK_EQUAL(CombinationGrid::enumerateLevelsWithSumWithoutBoundary(3, 7).size(), 15);

  sgpp::base::SBsplineBase basis1d;
  HeterogeneousBasis basis(3, basis1d);

  for (bool hasBoundary : {true, false}) {
    CombinationGrid combinationGrid = CombinationGrid::regular(3, 5, basis, hasBoundary);
    BOOST_CHECK_EQUAL(combinationGrid.getFullGrids().size(),
        (hasBoundary ? (21 + 15 + 10) : (15 + 10 + 6)));
    sgpp::base::GridStorage gridStorage(3);
    combinationGrid.combinePoints(gridStorage);
    BOOST_CHECK_EQUAL(gridStorage.getSize(), (hasBoundary ? 705 : 351));
  }

  /*std::cout << "BOOM1\n";
  CombinationGrid combinationGrid = CombinationGrid::regular(20, 6, basis, false);
  std::cout << "BOOM2: " << combinationGrid.getCoefficients().size() << "\n";
  sgpp::base::GridStorage gridStorage(20);
  combinationGrid.combinePoints(gridStorage);
  std::cout << "BOOM3: " << gridStorage.getSize() << "\n";*/

  /*for (const CombinationGrid& combinationGrid  : {
      CombinationGrid::regular(3, 5, basis), CombinationGrid::regular(3, 5, basis, false)}) {
    for (size_t i = 0; i < combinationGrid.getCoefficients().size(); i++) {
      const LevelVector& level = combinationGrid.getFullGrids()[i].getLevel();
      const double coefficient = combinationGrid.getCoefficients()[i];
      std::cout << level.size() << ": " << level[0] << ", " << level[1] << ", " << level[2] << " - " << coefficient << "\n";
    }

    std::cout << "\n";
  }*/
}

BOOST_AUTO_TEST_CASE(testOperationEvalCombinationGrid) {
  sgpp::base::SLinearBase basis1d;
  HeterogeneousBasis basis(2, basis1d);
  CombinationGrid combinationGrid = CombinationGrid::regular(2, 2, basis, false);
  OperationEvalCombinationGrid op(combinationGrid);
  std::vector<DataVector> surpluses = {
      DataVector{1.0, -2.0, -1.0}, DataVector{-0.5, 0.5, 1.0}, DataVector{5.0}};
  DataVector point{0.625, 0.375};
  // 1 * (1.0*(0*0.75) + -2.0*(0.5*0.75) + -1.0*(0.5*0.75))
  // + 1 * (-0.5*(0.75*0.5) + 0.5*(0.75*0.5) + 1.0*(0.75*0))
  // + (-1) * 5.0*0.75*0.75
  // = (0 - 0.75 - 0.375) + (-0.1875 + 0.1875 + 0) - 2.8125 = -1.125 + 0 - 2.8125 = -3.9375
  BOOST_CHECK_EQUAL(op.eval(surpluses, point), -3.9375);
}

BOOST_AUTO_TEST_CASE(testOperationUPFullGridLinear) {
  sgpp::base::SLinearBase basis1d;
  HeterogeneousBasis basis(2, basis1d);
  FullGrid fullGrid({2, 1}, basis);
  OperationPoleHierarchisationLinear operationPole;
  OperationUPFullGrid operation(fullGrid, operationPole);
  DataVector values{-0.5, 3.0, 0.25, 0.5, -1.0, 1.0, 5.0, 2.5,
                    -1.5, 0.0, 2.0, -1.0, 1.0, -2.0, -1.0};
  operation.apply(values);
  DataVector correctSurpluses{-0.5, 3.125, 1.0, 0.875, -1.0, 0.25, 2.9375, 1.25,
                              -2.1875, 1.0, 2.0, -2.5, 0.5, -2.0, -1.0};

  for (size_t i = 0; i < values.size(); i++) {
    BOOST_CHECK_CLOSE(values[i], correctSurpluses[i], 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(testOperationUPFullGridGeneral) {
  sgpp::base::SBsplineBase basis1d(3);
  HeterogeneousBasis basis(2, basis1d);
  FullGrid fullGrid({2, 1}, basis);
  std::vector<std::unique_ptr<OperationPole>> operationPole;
  OperationPoleHierarchisationGeneral::fromHeterogenerousBasis(basis, operationPole);
  OperationUPFullGrid operation(fullGrid, operationPole);
  DataVector values{-0.5, 3.0, 0.25, 0.5, -1.0, 1.0, 5.0, 2.5,
                    -1.5, 0.0, 2.0, -1.0, 1.0, -2.0, -1.0};
  operation.apply(values);
  DataVector correctSurpluses{-4.34548191029053, 11.29388598889826, -0.45638469432632633,
      8.30837174457185, -3.402942074462426, -4.7579856595053664, 10.752772401474628,
      4.220642134584388, -8.117821514507945, 5.795370207225638, 8.991769884961288,
      -15.38082204309674, 3.036651517372788, -7.852574819533291, -3.5702774351738715};

  for (size_t i = 0; i < values.size(); i++) {
    BOOST_CHECK_CLOSE(values[i], correctSurpluses[i], 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(testOperationUPCombinationGrid) {
  sgpp::base::SBsplineBase basis1d;
  HeterogeneousBasis basis(2, basis1d);
  CombinationGrid combinationGrid = CombinationGrid::regular(2, 1, basis);
  OperationPoleNodalisationBspline operationPole(3);
  OperationUPCombinationGrid operation(combinationGrid, operationPole);
  std::vector<DataVector> values = {
      DataVector{-0.5, 3.0, 0.25, 0.5, -1.0, 1.0}, DataVector{5.0, 2.5, -1.5, 0.0, 2.0, -1.0},
      DataVector{1.0, -2.0, -1.0, 0.5}};
  operation.apply(values);
  BOOST_CHECK_CLOSE(values[0][0], -3.83571429, 1e-6);
  BOOST_CHECK_CLOSE(values[0][1], 9.34285714, 1e-6);
  BOOST_CHECK_CLOSE(values[0][2], -2.33571429, 1e-6);
  BOOST_CHECK_CLOSE(values[0][3], 2.96785714, 1e-6);
  BOOST_CHECK_CLOSE(values[0][4], -5.87142857, 1e-6);
  BOOST_CHECK_CLOSE(values[0][5], 3.71785714, 1e-6);
  BOOST_CHECK_CLOSE(values[1][0], 12.66428571, 1e-6);
  BOOST_CHECK_CLOSE(values[1][1], 2.7, 1e-6);
  BOOST_CHECK_CLOSE(values[1][2], -8.65714286, 1e-6);
  BOOST_CHECK_CLOSE(values[1][3], 1.2, 1e-6);
  BOOST_CHECK_CLOSE(values[1][4], 7.56428571, 1e-6);
  BOOST_CHECK_CLOSE(values[1][5], -3.9, 1e-6);
  BOOST_CHECK_CLOSE(values[2][0], 4.56, 1e-6);
  BOOST_CHECK_CLOSE(values[2][1], -6.24, 1e-6);
  BOOST_CHECK_CLOSE(values[2][2], -3.84, 1e-6);
  BOOST_CHECK_CLOSE(values[2][3], 3.36, 1e-6);
}
