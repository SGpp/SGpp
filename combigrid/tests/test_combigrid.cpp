// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppCombigridModule

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/tools/Printer.hpp>

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/AdaptiveCombinationGridGenerator.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationEvalCombinationGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>
#include <sgpp/combigrid/operation/OperationPoleHierarchisationGeneral.hpp>
#include <sgpp/combigrid/operation/OperationPoleHierarchisationLinear.hpp>
#include <sgpp/combigrid/operation/OperationPoleNodalisationBspline.hpp>
#include <sgpp/combigrid/operation/OperationPoleNodalisationLinear.hpp>
#include <sgpp/combigrid/operation/OperationUPCombinationGrid.hpp>
#include <sgpp/combigrid/operation/OperationUPFullGrid.hpp>
#include <sgpp/combigrid/tools/IndexVectorRange.hpp>
#include <sgpp/combigrid/tools/LevelVectorTools.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <functional>
#include <memory>
#include <numeric>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::combigrid::AdaptiveCombinationGridGenerator;
using sgpp::combigrid::CombinationGrid;
using sgpp::combigrid::FullGrid;
using sgpp::combigrid::HeterogeneousBasis;
using sgpp::combigrid::IndexVector;
using sgpp::combigrid::IndexVectorRange;
using sgpp::combigrid::LevelVector;
using sgpp::combigrid::LevelVectorTools;
using sgpp::combigrid::OperationEvalCombinationGrid;
using sgpp::combigrid::OperationPole;
using sgpp::combigrid::OperationPoleHierarchisationGeneral;
using sgpp::combigrid::OperationPoleHierarchisationLinear;
using sgpp::combigrid::OperationPoleNodalisationBspline;
using sgpp::combigrid::OperationPoleNodalisationLinear;
using sgpp::combigrid::OperationUPCombinationGrid;
using sgpp::combigrid::OperationUPFullGrid;

// fix for clang (from https://stackoverflow.com/a/33755176)
#ifdef __clang__
#include <string>

namespace boost {
namespace unit_test {
namespace ut_detail {

std::string normalize_test_case_name(const_string name) {
  return ((name[0] == '&') ? std::string(name.begin() + 1, name.size() - 1)
                           : std::string(name.begin(), name.size()));
}

}  // namespace ut_detail
}  // namespace unit_test
}  // namespace boost
#endif

BOOST_AUTO_TEST_CASE(testFullGrid) {
  BOOST_CHECK_EQUAL(FullGrid().getDimension(), 0);

  sgpp::base::SBsplineBase basis1d;
  const HeterogeneousBasis basis(3, basis1d);
  const FullGrid fullGrid({2, 3, 2}, basis);

  BOOST_CHECK(fullGrid == FullGrid({2, 3, 2}, basis));
  BOOST_CHECK(!(fullGrid == FullGrid({2, 2, 2}, basis)));

  BOOST_CHECK(fullGrid != FullGrid({2, 2, 2}, basis));
  BOOST_CHECK(!(fullGrid != FullGrid({2, 3, 2}, basis)));

  BOOST_CHECK_EQUAL(fullGrid.getMinIndex(0), 0);
  BOOST_CHECK_EQUAL(fullGrid.getMinIndex(1), 0);
  BOOST_CHECK_EQUAL(fullGrid.getMinIndex(2), 0);

  BOOST_CHECK_EQUAL(fullGrid.getMaxIndex(0), 4);
  BOOST_CHECK_EQUAL(fullGrid.getMaxIndex(1), 8);
  BOOST_CHECK_EQUAL(fullGrid.getMaxIndex(2), 4);

  BOOST_CHECK_EQUAL(fullGrid.getNumberOfIndexVectors(0), 5);
  BOOST_CHECK_EQUAL(fullGrid.getNumberOfIndexVectors(1), 9);
  BOOST_CHECK_EQUAL(fullGrid.getNumberOfIndexVectors(2), 5);
  BOOST_CHECK_EQUAL(fullGrid.getNumberOfIndexVectors(), 225);
}

BOOST_AUTO_TEST_CASE(testHeterogeneousBasis) {
  sgpp::base::SBsplineBase basis1d1(1);
  sgpp::base::SBsplineBase basis1d3(3);
  sgpp::base::SBsplineBase basis1d5(5);

  BOOST_CHECK(HeterogeneousBasis(3, basis1d3) ==
              HeterogeneousBasis({&basis1d3, &basis1d3, &basis1d3}));

  {
    sgpp::combigrid::level_t level;
    sgpp::combigrid::index_t index;

    level = 6;
    index = 0;
    HeterogeneousBasis::hierarchizeLevelIndex(level, index);
    BOOST_CHECK_EQUAL(level, 0);
    BOOST_CHECK_EQUAL(index, 0);

    level = 4;
    index = 16;
    HeterogeneousBasis::hierarchizeLevelIndex(level, index);
    BOOST_CHECK_EQUAL(level, 0);
    BOOST_CHECK_EQUAL(index, 1);

    level = 4;
    index = 12;
    HeterogeneousBasis::hierarchizeLevelIndex(level, index);
    BOOST_CHECK_EQUAL(level, 2);
    BOOST_CHECK_EQUAL(index, 3);

    level = 6;
    index = 11;
    HeterogeneousBasis::hierarchizeLevelIndex(level, index);
    BOOST_CHECK_EQUAL(level, 6);
    BOOST_CHECK_EQUAL(index, 11);
  }

  const HeterogeneousBasis basisHierarchical({&basis1d1, &basis1d3, &basis1d5}, true);
  const HeterogeneousBasis basisNodal({&basis1d1, &basis1d3, &basis1d5}, false);

  // 0.8 * 0.53866667 * 0.53037333 = 0.2285555
  BOOST_CHECK_CLOSE(basisHierarchical.eval({1, 2, 3}, {1, 2, 6}, DataVector({0.6, 0.7, 0.8})),
                    0.2285555, 1e-4);
  // 0.8 * 0.28266667 * 0.47554667 = 0.1075370
  BOOST_CHECK_CLOSE(basisNodal.eval({1, 2, 3}, {1, 2, 6}, DataVector({0.6, 0.7, 0.8})), 0.1075370,
                    1e-4);
}

BOOST_AUTO_TEST_CASE(testCombinationGrid) {
  BOOST_CHECK_EQUAL(LevelVectorTools::generateDiagonalWithBoundary(3, 5).size(), 21);
  BOOST_CHECK_EQUAL(LevelVectorTools::generateDiagonalWithoutBoundary(3, 7).size(), 15);

  sgpp::base::SBsplineBase basis1d;
  HeterogeneousBasis basis(3, basis1d);

  for (bool hasBoundary : {true, false}) {
    std::vector<LevelVector> subspaceLevels;

    if (hasBoundary) {
      subspaceLevels = {
          {0, 0, 5}, {0, 1, 4}, {0, 2, 3}, {0, 3, 2}, {0, 4, 1}, {0, 5, 0}, {1, 0, 4}, {1, 1, 3},
          {1, 2, 2}, {1, 3, 1}, {1, 4, 0}, {2, 0, 3}, {2, 1, 2}, {2, 2, 1}, {2, 3, 0}, {3, 0, 2},
          {3, 1, 1}, {3, 2, 0}, {4, 0, 1}, {4, 1, 0}, {5, 0, 0}, {0, 0, 4}, {0, 1, 3}, {0, 2, 2},
          {0, 3, 1}, {0, 4, 0}, {1, 0, 3}, {1, 1, 2}, {1, 2, 1}, {1, 3, 0}, {2, 0, 2}, {2, 1, 1},
          {2, 2, 0}, {3, 0, 1}, {3, 1, 0}, {4, 0, 0}, {0, 0, 3}, {0, 1, 2}, {0, 2, 1}, {0, 3, 0},
          {1, 0, 2}, {1, 1, 1}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}, {3, 0, 0}, {0, 0, 2}, {0, 1, 1},
          {0, 2, 0}, {1, 0, 1}, {1, 1, 0}, {2, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 0}};
    } else {
      subspaceLevels = {{1, 1, 5}, {1, 2, 4}, {1, 3, 3}, {1, 4, 2}, {1, 5, 1}, {2, 1, 4},
                        {2, 2, 3}, {2, 3, 2}, {2, 4, 1}, {3, 1, 3}, {3, 2, 2}, {3, 3, 1},
                        {4, 1, 2}, {4, 2, 1}, {5, 1, 1}, {1, 1, 4}, {1, 2, 3}, {1, 3, 2},
                        {1, 4, 1}, {2, 1, 3}, {2, 2, 2}, {2, 3, 1}, {3, 1, 2}, {3, 2, 1},
                        {4, 1, 1}, {1, 1, 3}, {1, 2, 2}, {1, 3, 1}, {2, 1, 2}, {2, 2, 1},
                        {3, 1, 1}, {1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {1, 1, 1}};
    }

    std::vector<CombinationGrid> combinationGrids = {
        CombinationGrid::fromRegularSparse(3, 5, basis, hasBoundary),
        CombinationGrid::fromSubspaces(subspaceLevels, basis, hasBoundary)};

    for (const CombinationGrid& combinationGrid : combinationGrids) {
      BOOST_CHECK_EQUAL(combinationGrid.getFullGrids().size(),
                        (hasBoundary ? (21 + 15 + 10) : (15 + 10 + 6)));
      sgpp::base::GridStorage gridStorage(3);
      combinationGrid.combinePoints(gridStorage);
      BOOST_CHECK_EQUAL(gridStorage.getSize(), (hasBoundary ? 705 : 351));
    }

    const std::vector<FullGrid>& fullGrids0 = combinationGrids[0].getFullGrids();
    const std::vector<FullGrid>& fullGrids1 = combinationGrids[1].getFullGrids();
    const DataVector& coefficients0 = combinationGrids[0].getCoefficients();
    const DataVector& coefficients1 = combinationGrids[1].getCoefficients();

    for (size_t i = 0; i < fullGrids1.size(); i++) {
      const FullGrid& fullGrid = fullGrids1[i];
      auto it = std::find(fullGrids0.begin(), fullGrids0.end(), fullGrid);
      BOOST_CHECK(it != fullGrids0.end());
      BOOST_CHECK_EQUAL(coefficients0[it - fullGrids0.begin()], coefficients1[i]);
    }
  }

  CombinationGrid combinationGrid = CombinationGrid::fromRegularSparse(2, 2, basis);
  BOOST_CHECK_EQUAL(combinationGrid.combineValues(DataVector({1.0, 0.5, -0.5, -2.0, 4.0})), -1.0);

  sgpp::base::GridStorage gridStorage(2);
  combinationGrid.combinePoints(gridStorage);
  BOOST_CHECK_EQUAL(gridStorage.getSize(), 17);
  const std::vector<DataVector> values = {
      // 31013
      // |   |
      // |   |
      // |   |
      // 23221
      DataVector({2.0, 3.0, 2.0, 2.0, 1.0, 3.0, 1.0, 0.0, 1.0, 3.0}),
      // 2-2-1
      // |   |
      // 0 2 1
      // |   |
      // 2-1-0
      DataVector({2.0, 1.0, 0.0, 0.0, 2.0, 1.0, 2.0, 2.0, 1.0}),
      // 1---0
      // 3   0
      // 5   3
      // 0   1
      // 1---2
      DataVector({1.0, 2.0, 0.0, 1.0, 5.0, 3.0, 3.0, 0.0, 1.0, 0.0}),
      // 0-1-1
      // |   |
      // |   |
      // |   |
      // 3-1-1
      DataVector({3.0, 1.0, 1.0, 0.0, 1.0, 1.0}),
      // 2---0
      // |   |
      // 2   1
      // |   |
      // 2---1
      DataVector({2.0, 1.0, 2.0, 1.0, 2.0, 0.0})};
  // 41113
  // 3   0
  // 3 2 3
  // 0   1
  // 03221
  const DataVector correctResult(
      {0.0, 3.0, 2.0, 2.0, 1.0, 0.0, 1.0, 3.0, 2.0, 3.0, 3.0, 0.0, 4.0, 1.0, 1.0, 1.0, 3.0});
  const std::vector<DataVector> correctDistributedValues = {
      // 41113
      // |   |
      // |   |
      // |   |
      // 03221
      DataVector({0.0, 3.0, 2.0, 2.0, 1.0, 4.0, 1.0, 1.0, 1.0, 3.0}),
      // 4-1-3
      // |   |
      // 3 2 3
      // |   |
      // 0-2-1
      DataVector({0.0, 2.0, 1.0, 3.0, 2.0, 3.0, 4.0, 1.0, 3.0}),
      // 4---3
      // 3   0
      // 3   3
      // 0   1
      // 0---1
      DataVector({0.0, 1.0, 0.0, 1.0, 3.0, 3.0, 3.0, 0.0, 4.0, 3.0}),
      // 4-1-3
      // |   |
      // |   |
      // |   |
      // 0-2-1
      DataVector({0.0, 2.0, 1.0, 4.0, 1.0, 3.0}),
      // 4---3
      // |   |
      // 3   3
      // |   |
      // 0---1
      DataVector({0.0, 1.0, 3.0, 3.0, 4.0, 3.0})};

  DataVector result;
  combinationGrid.combineSparseGridValues(gridStorage, values, result);

  std::vector<size_t> order(gridStorage.getSize());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&gridStorage](size_t a, size_t b) {
    const double a1 = gridStorage.getPoint(a).getStandardCoordinate(1);
    const double b1 = gridStorage.getPoint(b).getStandardCoordinate(1);

    if (a1 != b1) {
      return (a1 < b1);
    } else {
      const double a0 = gridStorage.getPoint(a).getStandardCoordinate(0);
      const double b0 = gridStorage.getPoint(b).getStandardCoordinate(0);
      return (a0 < b0);
    }
  });

  for (size_t k = 0; k < gridStorage.getSize(); k++) {
    BOOST_CHECK_EQUAL(result[order[k]], correctResult[k]);
  }

  std::vector<DataVector> distributedValues;
  combinationGrid.distributeValuesToFullGrids(gridStorage, result, distributedValues);
  BOOST_CHECK_EQUAL(distributedValues.size(), correctDistributedValues.size());

  for (size_t i = 0; i < distributedValues.size(); i++) {
    BOOST_CHECK_EQUAL_COLLECTIONS(distributedValues[i].begin(), distributedValues[i].end(),
                                  correctDistributedValues[i].begin(),
                                  correctDistributedValues[i].end());
  }
}

BOOST_AUTO_TEST_CASE(testIndexVectorRangeIterator) {
  sgpp::base::SBsplineBase basis1d;
  const HeterogeneousBasis basis(2, basis1d);
  const FullGrid fullGrid({1, 2}, basis);

  for (size_t j = 0; j < 2; j++) {
    IndexVectorRange range;
    std::vector<IndexVector> correctIndices;

    if (j == 0) {
      range = IndexVectorRange(fullGrid);
      correctIndices = {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}, {0, 2}, {1, 2},
                        {2, 2}, {0, 3}, {1, 3}, {2, 3}, {0, 4}, {1, 4}, {2, 4}};
      BOOST_CHECK_EQUAL(range.find({2, 1}), 5);
    } else {
      range = IndexVectorRange({3, 6, 1}, {6, 8, 2});
      correctIndices = {{3, 6, 1}, {4, 6, 1}, {5, 6, 1}, {6, 6, 1}, {3, 7, 1}, {4, 7, 1},
                        {5, 7, 1}, {6, 7, 1}, {3, 8, 1}, {4, 8, 1}, {5, 8, 1}, {6, 8, 1},
                        {3, 6, 2}, {4, 6, 2}, {5, 6, 2}, {6, 6, 2}, {3, 7, 2}, {4, 7, 2},
                        {5, 7, 2}, {6, 7, 2}, {3, 8, 2}, {4, 8, 2}, {5, 8, 2}, {6, 8, 2}};
      BOOST_CHECK_EQUAL(range.find({5, 7, 2}), 18);
    }

    BOOST_CHECK_EQUAL(range.end() - range.begin(), correctIndices.size());
    size_t i = 0;

    for (const IndexVector& index : range) {
      BOOST_CHECK_EQUAL_COLLECTIONS(index.begin(), index.end(), correctIndices[i].begin(),
                                    correctIndices[i].end());
      i++;
    }
  }

  DataMatrix points;
  IndexVectorRange::getPoints(fullGrid, points);
  const DataMatrix correctPoints(
      {0.0, 0.0, 0.5, 0.0, 1.0,  0.0, 0.0,  0.25, 0.5,  0.25, 1.0, 0.25, 0.0, 0.5, 0.5,
       0.5, 1.0, 0.5, 0.0, 0.75, 0.5, 0.75, 1.0,  0.75, 0.0,  1.0, 0.5,  1.0, 1.0, 1.0},
      15);
  BOOST_CHECK_EQUAL_COLLECTIONS(points.begin(), points.end(), correctPoints.begin(),
                                correctPoints.end());
}

BOOST_AUTO_TEST_CASE(testOperationEvalCombinationGrid) {
  sgpp::base::SLinearBase basis1d;
  const HeterogeneousBasis basis(2, basis1d, false);
  const CombinationGrid combinationGrid = CombinationGrid::fromRegularSparse(2, 2, basis, false);
  OperationEvalCombinationGrid op(combinationGrid);
  const std::vector<DataVector> surpluses = {DataVector{1.0, -2.0, -1.0},
                                             DataVector{-0.5, 0.5, 1.0}, DataVector{5.0}};
  DataVector point({0.625, 0.375});
  // 1 * (1.0*(0*0.75) + -2.0*(0.5*0.75) + -1.0*(0.5*0.75))
  // + 1 * (-0.5*(0.75*0.5) + 0.5*(0.75*0.5) + 1.0*(0.75*0))
  // + (-1) * 5.0*0.75*0.75
  // = (0 - 0.75 - 0.375) + (-0.1875 + 0.1875 + 0) - 2.8125 = -1.125 + 0 - 2.8125 = -3.9375
  BOOST_CHECK_EQUAL(op.eval(surpluses, point), -3.9375);

  DataMatrix points({0.625, 0.375, 0.625, 0.375}, 2);
  DataVector result;
  op.multiEval(surpluses, points, result);
  BOOST_CHECK_EQUAL(result[0], -3.9375);
  BOOST_CHECK_EQUAL(result[1], -3.9375);
}

BOOST_AUTO_TEST_CASE(testOperationUPFullGridLinear) {
  sgpp::base::SLinearBase basis1d;
  const HeterogeneousBasis basis(2, basis1d);
  const FullGrid fullGrid({2, 1}, basis);
  OperationPoleHierarchisationLinear operationPole;
  OperationUPFullGrid operation(fullGrid, operationPole);
  DataVector values{-0.5, 3.0, 0.25, 0.5,  -1.0, 1.0,  5.0, 2.5,
                    -1.5, 0.0, 2.0,  -1.0, 1.0,  -2.0, -1.0};
  operation.apply(values);
  const DataVector correctSurpluses{-0.5,    3.125, 1.0, 0.875, -1.0, 0.25, 2.9375, 1.25,
                                    -2.1875, 1.0,   2.0, -2.5,  0.5,  -2.0, -1.0};

  for (size_t i = 0; i < values.size(); i++) {
    BOOST_CHECK_CLOSE(values[i], correctSurpluses[i], 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(testOperationUPFullGridGeneral) {
  for (bool isBasisHierarchical : {true, false}) {
    sgpp::base::SBsplineBase basis1d(3);
    const HeterogeneousBasis basis(2, basis1d, isBasisHierarchical);
    const FullGrid fullGrid({2, 1}, basis);
    std::vector<std::unique_ptr<OperationPole>> operationPole;
    OperationPoleHierarchisationGeneral::fromHeterogenerousBasis(basis, operationPole);
    OperationUPFullGrid operation(fullGrid, operationPole);
    DataVector values{-0.5, 3.0, 0.25, 0.5,  -1.0, 1.0,  5.0, 2.5,
                      -1.5, 0.0, 2.0,  -1.0, 1.0,  -2.0, -1.0};
    operation.apply(values);
    DataVector correctSurpluses;

    if (isBasisHierarchical) {
      correctSurpluses = DataVector{-4.34548191029053,  11.29388598889826,  -0.45638469432632633,
                                    8.30837174457185,   -3.402942074462426, -4.7579856595053664,
                                    10.752772401474628, 4.220642134584388,  -8.117821514507945,
                                    5.795370207225638,  8.991769884961288,  -15.38082204309674,
                                    3.036651517372788,  -7.852574819533291, -3.5702774351738715};
    } else {
      correctSurpluses = DataVector{-2.860096153846154, 5.333241758241758,   -3.044299450549451,
                                    3.4689560439560445, -3.4386675824175827, -1.101923076923078,
                                    10.83626373626374,  4.042582417582418,   -4.506593406593407,
                                    2.4123626373626372, 5.998557692307693,   -7.601373626373627,
                                    3.8355082417582422, -4.365659340659341,  -1.4800137362637362};
    }

    for (size_t i = 0; i < values.size(); i++) {
      BOOST_CHECK_CLOSE(values[i], correctSurpluses[i], 1e-8);
    }
  }
}

BOOST_AUTO_TEST_CASE(testOperationUPCombinationGrid) {
  sgpp::base::SBsplineBase basis1d;
  const HeterogeneousBasis basis(2, basis1d);
  const CombinationGrid combinationGrid = CombinationGrid::fromRegularSparse(2, 1, basis);
  const std::vector<DataVector> originalValues = {DataVector{-0.5, 3.0, 0.25, 0.5, -1.0, 1.0},
                                                  DataVector{5.0, 2.5, -1.5, 0.0, 2.0, -1.0},
                                                  DataVector{1.0, -2.0, -1.0, 0.5}};

  {
    OperationPoleNodalisationBspline operationPole(3);
    OperationUPCombinationGrid operation(combinationGrid, operationPole);
    std::vector<DataVector> values = originalValues;

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

  {
    OperationPoleNodalisationLinear operationPole;
    OperationUPCombinationGrid operation(combinationGrid, operationPole);
    std::vector<DataVector> values = originalValues;

    operation.apply(values);

    BOOST_CHECK_EQUAL(values[0][0], -0.5);
    BOOST_CHECK_EQUAL(values[0][1], 3.0);
    BOOST_CHECK_EQUAL(values[0][2], 0.25);
    BOOST_CHECK_EQUAL(values[0][3], 0.5);
    BOOST_CHECK_EQUAL(values[0][4], -1.0);
    BOOST_CHECK_EQUAL(values[0][5], 1.0);
    BOOST_CHECK_EQUAL(values[1][0], 5.0);
    BOOST_CHECK_EQUAL(values[1][1], 2.5);
    BOOST_CHECK_EQUAL(values[1][2], -1.5);
    BOOST_CHECK_EQUAL(values[1][3], 0.0);
    BOOST_CHECK_EQUAL(values[1][4], 2.0);
    BOOST_CHECK_EQUAL(values[1][5], -1.0);
    BOOST_CHECK_EQUAL(values[2][0], 1.0);
    BOOST_CHECK_EQUAL(values[2][1], -2.0);
    BOOST_CHECK_EQUAL(values[2][2], -1.0);
    BOOST_CHECK_EQUAL(values[2][3], 0.5);
  }
}

namespace std {
// needed for BOOST_CHECK_EQUAL_COLLECTIONS in next test
using sgpp::base::operator<<;
}  // namespace std

BOOST_AUTO_TEST_CASE(testMakeDownwardClosed) {
  std::vector<LevelVector> subspaces = {LevelVector{0, 0, 1}, LevelVector{0, 2, 1},
                                        LevelVector{1, 0, 3}};

  std::vector<LevelVector> downwardClosedSetSolution = {
      LevelVector{0, 0, 0}, LevelVector{1, 0, 0}, LevelVector{0, 1, 0}, LevelVector{0, 2, 0},
      LevelVector{0, 0, 1}, LevelVector{1, 0, 1}, LevelVector{0, 1, 1}, LevelVector{0, 2, 1},
      LevelVector{0, 0, 2}, LevelVector{1, 0, 2}, LevelVector{0, 0, 3}, LevelVector{1, 0, 3},
  };

  auto downwardClosedSet = LevelVectorTools::makeDownwardClosed(LevelVector{0, 0, 0}, subspaces);

  BOOST_CHECK_EQUAL_COLLECTIONS(downwardClosedSetSolution.begin(), downwardClosedSetSolution.end(),
                                downwardClosedSet.begin(), downwardClosedSet.end());
}

BOOST_AUTO_TEST_CASE(testAdaptiveCombinationGridGenerator) {
  using sgpp::base::operator<<;
  for (bool hasBoundary : {true, false}) {
    sgpp::base::SBsplineBase basis1d;
    HeterogeneousBasis basis(3, basis1d);

    auto combinationGrid = CombinationGrid::fromRegularSparse(3, 5, basis, hasBoundary);

    auto adaptiveCombinationGridGenerator = AdaptiveCombinationGridGenerator::fromCombinationGrid(
        combinationGrid, std::plus<double>(),
        std::unique_ptr<sgpp::combigrid::RelevanceCalculator>(
            new sgpp::combigrid::WeightedRelevanceCalculator()),
        std::unique_ptr<sgpp::combigrid::PriorityEstimator>(
            new sgpp::combigrid::AveragingPriorityEstimator()));
    if (hasBoundary) {
      BOOST_CHECK_EQUAL(adaptiveCombinationGridGenerator.getMinimumLevelVector(),
                        LevelVector(3, 0));
    } else {
      BOOST_CHECK_EQUAL(adaptiveCombinationGridGenerator.getMinimumLevelVector(),
                        LevelVector(3, 1));
    }

    // feed "known values" for the initial combinationGrid to the adaptiveCombinationGridGenerator
    std::vector<LevelVector> subspaces{};
    subspaces.resize(combinationGrid.getFullGrids().size());
    std::transform(combinationGrid.getFullGrids().begin(), combinationGrid.getFullGrids().end(),
                   subspaces.begin(),
                   [](const FullGrid& fg) -> LevelVector { return fg.getLevel(); });
    std::sort(subspaces.begin(), subspaces.end());
    for (const auto& subspace : subspaces) {
      adaptiveCombinationGridGenerator.setQoIInformation(subspace, 1.);
      bool notAdapted = !adaptiveCombinationGridGenerator.adaptNextLevelVector();
      // those should all already be in the old set, from the way the generator was constructed
      BOOST_CHECK(notAdapted);
    }
    bool notAdapted = !adaptiveCombinationGridGenerator.adaptAllKnown();
    BOOST_CHECK(notAdapted);

    // adapt three times, should be regular
    for (size_t i = 0; i < 3; i++) {
      auto activeSet = adaptiveCombinationGridGenerator.getActiveSet();
      for (auto& active : activeSet) {
        adaptiveCombinationGridGenerator.setQoIInformation(active, 1.);
      }
      BOOST_CHECK(adaptiveCombinationGridGenerator.adaptAllKnown());
    }

    // assert that we get the same subspaces we would get when taking a larger regular
    // combinationGrid
    auto adaptedCombinationGrid = adaptiveCombinationGridGenerator.getCombinationGrid(basis);

    std::vector<LevelVector> adaptedSubspaces{};
    adaptedSubspaces.resize(adaptedCombinationGrid.getFullGrids().size());
    std::transform(adaptedCombinationGrid.getFullGrids().begin(),
                   adaptedCombinationGrid.getFullGrids().end(), adaptedSubspaces.begin(),
                   [](const FullGrid& fg) -> LevelVector { return fg.getLevel(); });
    std::sort(adaptedSubspaces.begin(), adaptedSubspaces.end());

    auto largerCombinationGrid = CombinationGrid::fromRegularSparse(3, 8, basis, hasBoundary);
    std::vector<LevelVector> largerSubspaces{};
    largerSubspaces.resize(largerCombinationGrid.getFullGrids().size());
    std::transform(largerCombinationGrid.getFullGrids().begin(),
                   largerCombinationGrid.getFullGrids().end(), largerSubspaces.begin(),
                   [](const FullGrid& fg) -> LevelVector { return fg.getLevel(); });
    std::sort(largerSubspaces.begin(), largerSubspaces.end());

    BOOST_CHECK_EQUAL_COLLECTIONS(adaptedSubspaces.begin(), adaptedSubspaces.end(),
                                  largerSubspaces.begin(), largerSubspaces.end());
  }
}
