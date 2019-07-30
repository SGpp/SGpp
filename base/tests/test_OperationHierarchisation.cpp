// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>

using sgpp::base::DataVector;
using sgpp::base::BoundingBox1D;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::OperationEval;
using sgpp::base::OperationHierarchisation;
using sgpp::base::Stretching;
using sgpp::base::Stretching1D;

void testHierarchisationDehierarchisation(sgpp::base::Grid& grid, size_t level,
                                          double (*func)(DataVector&), double tolerance = 1e-12,
                                          bool naiveOp = false, bool arbitraryBoundary = false) {
  grid.getGenerator().regular(level);
  GridStorage& gridStore = grid.getStorage();
  size_t dim = gridStore.getDimension();

  DataVector node_values = DataVector(gridStore.getSize());
  DataVector coords = DataVector(dim);

  for (size_t n = 0; n < gridStore.getSize(); n++) {
    gridStore.getCoordinates(gridStore[n], coords);
    node_values[n] = func(coords);
  }

  DataVector alpha = DataVector(node_values);
  std::unique_ptr<OperationHierarchisation> hierarchisation;
  if (arbitraryBoundary) {
    hierarchisation.reset(sgpp::op_factory::createOperationArbitraryBoundaryHierarchisation(grid));
  } else {
    hierarchisation.reset(sgpp::op_factory::createOperationHierarchisation(grid));
  }
  hierarchisation->doHierarchisation(alpha);

  std::unique_ptr<OperationEval> op;

  if (naiveOp == true) {
    op.reset(sgpp::op_factory::createOperationEvalNaive(grid));
  } else {
    op.reset(sgpp::op_factory::createOperationEval(grid));
  }

  for (size_t n = 0; n < gridStore.getSize(); n++) {
    gridStore.getCoordinates(gridStore[n], coords);
    double eval = op->eval(alpha, coords);

    if (std::abs(node_values[n]) < 1e-10) {
      BOOST_CHECK_SMALL(eval - node_values[n], tolerance);
    } else {
      BOOST_CHECK_CLOSE(eval, node_values[n], tolerance);
    }
  }

  DataVector node_values_back = DataVector(alpha);
  hierarchisation->doDehierarchisation(node_values_back);

  for (size_t n = 0; n < gridStore.getSize(); n++) {
    if (std::abs(node_values[n]) < 1e-10) {
      BOOST_CHECK_SMALL(node_values_back[n] - node_values[n], tolerance);
    } else {
      BOOST_CHECK_CLOSE(node_values_back[n], node_values[n], tolerance);
    }
  }
}

double parabola(DataVector& input) {
  double result = 1.;

  for (size_t i = 0; i < input.getSize(); i++) {
    result *= input[i] * (1. - input[i]) * 4.;
  }

  return result;
}

double parabolaBoundary(DataVector& input) {
  double result = 1.;

  for (size_t i = 0; i < input.getSize(); i++) {
    result *= 0.25 * (input[i] - 0.7) * (input[i] - 0.7) + 2;
  }

  return result;
}

BOOST_AUTO_TEST_SUITE(testHierarchization)

BOOST_AUTO_TEST_CASE(testHierarchisationLinear) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModLinear) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationLinearClenshawCurtis) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createLinearClenshawCurtisGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-13, true);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModLinearClenshawCurtis) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createModLinearClenshawCurtisGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationLinearClenshawCurtisBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createLinearClenshawCurtisBoundaryGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-13, true);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModLinearWithBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationTruncatedBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim, 0));
    testHierarchisationDehierarchisation(*grid, level, &parabola);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationLdBoundary) {
  int level = 5;

  for (int boundaryLevel = 2; boundaryLevel < level; boundaryLevel++) {
    for (int dim = 1; dim < 4; dim++) {
      std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim, boundaryLevel));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPrewavelet) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    std::unique_ptr<Grid> grid(Grid::createPrewaveletGrid(dim));
    testHierarchisationDehierarchisation(*grid, level, &parabola);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationFundamentalNakSplineBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    for (int degree = 1; degree < 6; degree += 2) {
      std::unique_ptr<Grid> grid(Grid::createFundamentalNakSplineBoundaryGrid(dim, degree));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-9, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationFundamentalSpline) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    for (int degree = 1; degree < 6; degree += 2) {
      std::unique_ptr<Grid> grid(Grid::createFundamentalSplineGrid(dim, degree));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-9, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationFundamentalSplineBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    for (int degree = 1; degree < 6; degree += 2) {
      std::unique_ptr<Grid> grid(Grid::createFundamentalSplineBoundaryGrid(dim, degree));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-9, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModFundamentalSpline) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    for (int degree = 1; degree < 6; degree += 2) {
      std::unique_ptr<Grid> grid(Grid::createModFundamentalSplineGrid(dim, degree));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-9, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPoly) {
  int level = 5;
  int degree[4] = {2, 3, 5, 8};

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      std::unique_ptr<Grid> grid(Grid::createPolyGrid(dim, degree[i]));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPolyTruncatedBoundary) {
  int level = 5;
  int degree[4] = {2, 3, 5, 8};

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      std::unique_ptr<Grid> grid(Grid::createPolyBoundaryGrid(dim, degree[i]));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModPoly) {
  int level = 5;
  int degree[4] = {2, 3, 5, 8};

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      std::unique_ptr<Grid> grid(Grid::createModPolyGrid(dim, degree[i]));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-12, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPolyClenshawCurtisBoundary) {
  int level = 5;
  int degree[4] = {2, 3, 5, 8};

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      std::unique_ptr<Grid> grid(Grid::createPolyClenshawCurtisBoundaryGrid(dim, degree[i]));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-13, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPolyClenshawCurtis) {
  int level = 5;
  int degree[4] = {2, 3, 5, 8};

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      std::unique_ptr<Grid> grid(Grid::createPolyClenshawCurtisGrid(dim, degree[i]));
      testHierarchisationDehierarchisation(*grid, level, &parabola, 1e-13, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationStretchedTruncatedBoundary1D) {
  int dim = 1;
  int level = 5;
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 0.;
  str1d.xsi = 10;
  BoundingBox1D dimBound;
  dimBound.leftBoundary = 0.00001;
  dimBound.rightBoundary = 1.;
  Stretching stretch({dimBound}, {str1d});
  std::unique_ptr<Grid> grid(Grid::createLinearStretchedBoundaryGrid(dim));
  grid->getStorage().setStretching(stretch);
  testHierarchisationDehierarchisation(*grid, level, &parabolaBoundary, 1e-13, false);
}

BOOST_AUTO_TEST_CASE(testHierarchisationStretchedTruncatedBoundary3D) {
  int dim = 3;
  int level = 5;
  Stretching1D str1d;  // = new Stretching1D();
  str1d.type = "sinh";
  str1d.x_0 = 1.;
  str1d.xsi = 10;
  BoundingBox1D dimBound;  // = new BoundingBox1D();
  dimBound.leftBoundary = 0.001;
  dimBound.rightBoundary = 1.;

  std::vector<Stretching1D> stretch_vec;
  stretch_vec.push_back(str1d);
  stretch_vec.push_back(str1d);
  stretch_vec.push_back(str1d);
  std::vector<BoundingBox1D> dimBound_vec;
  dimBound_vec.push_back(dimBound);
  dimBound_vec.push_back(dimBound);
  dimBound_vec.push_back(dimBound);

  Stretching stretch(dimBound_vec, stretch_vec);
  std::unique_ptr<Grid> grid(Grid::createLinearStretchedBoundaryGrid(dim));
  grid->getStorage().setStretching(stretch);
  testHierarchisationDehierarchisation(*grid, level, &parabolaBoundary, 1e-12, false);
}

BOOST_AUTO_TEST_SUITE_END()
