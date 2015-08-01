#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

using namespace SGPP::base;

void testHierarchisationDehierarchisation(SGPP::base::Grid* grid, size_t level, SGPP::float_t (*func)(DataVector&), SGPP::float_t tolerance = 0.0, bool naiveOp = false, bool doStretch = false) {
  GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  GridStorage* gridStore = grid->getStorage();
  size_t dim = gridStore->dim();
  Stretching* stretch = NULL;

  if (doStretch)
    stretch = gridStore->getStretching();

  DataVector node_values = DataVector(gridStore->size());
  DataVector coords = DataVector(dim);

  for (size_t n = 0; n < gridStore->size(); n++) {
    if (doStretch) {
      gridStore->get(n)->getCoordsStretching(coords, *stretch);
    } else {
      gridStore->get(n)->getCoords(coords);
    }

    node_values[n] = func(coords);
  }

  DataVector alpha = DataVector(node_values);
  OperationHierarchisation* hierarchisation = SGPP::op_factory::createOperationHierarchisation(*grid);
  hierarchisation->doHierarchisation(alpha);

  if (naiveOp == true) {
    OperationNaiveEval* op = SGPP::op_factory::createOperationNaiveEval(*grid);

    for (size_t n = 0; n < gridStore->size(); n++) {
      if (doStretch) {
        gridStore->get(n)->getCoordsStretching(coords, *stretch);
      } else {
        gridStore->get(n)->getCoords(coords);
      }

      SGPP::float_t eval = op->eval(alpha, coords);
      BOOST_CHECK_CLOSE(eval, node_values[n], tolerance);
    }
  } else {
    OperationEval* op = SGPP::op_factory::createOperationEval(*grid);

    for (size_t n = 0; n < gridStore->size(); n++) {
      if (doStretch) {
        gridStore->get(n)->getCoordsStretching(coords, *stretch);
      } else {
        gridStore->get(n)->getCoords(coords);
      }

      SGPP::float_t eval = op->eval(alpha, coords);
      BOOST_CHECK_CLOSE(eval, node_values[n], tolerance);
    }
  }

  DataVector node_values_back = DataVector(alpha);
  hierarchisation->doDehierarchisation(node_values_back);

  for (size_t n = 0; n < gridStore->size(); n++) {
    BOOST_CHECK_CLOSE(node_values_back[n], node_values[n], tolerance);
  }
}



SGPP::float_t parabola(DataVector& input) {
  SGPP::float_t result = 1.;

  for (size_t i = 0; i < input.getSize(); i++) {
    result *= input[i] * (1. - input[i]) * 4.;
  }

  return result;
}

SGPP::float_t parabolaBoundary(DataVector& input) {
  SGPP::float_t result = 1.;

  for (size_t i = 0; i < input.getSize(); i++) {
    result *= 0.25 * (input[i] - 0.7) * (input[i] - 0.7) + 2;
  }

  return result;
}

BOOST_AUTO_TEST_SUITE(testHierarchization)

BOOST_AUTO_TEST_CASE(testHierarchisationLinear) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    Grid* grid = Grid::createLinearGrid(dim);
    testHierarchisationDehierarchisation(grid, level, &parabola);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModLinear) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    Grid* grid = Grid::createModLinearGrid(dim);
    testHierarchisationDehierarchisation(grid, level, &parabola);
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModLinearWithBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    Grid* grid = Grid::createModLinearGrid(dim);
#if USE_DOUBLE_PRECISION
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-12, true);
#else
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-4, true);
#endif
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationTruncatedBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    Grid* grid = Grid::createLinearTruncatedBoundaryGrid(dim);
#if USE_DOUBLE_PRECISION
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-12, true);
#else
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-4, true);
#endif
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationBoundary) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    Grid* grid = Grid::createLinearBoundaryGrid(dim);
#if USE_DOUBLE_PRECISION
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-12 );
#else
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-4 );
#endif
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPrewavelet) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    Grid* grid = Grid::createPrewaveletGrid(dim);
#if USE_DOUBLE_PRECISION
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-12 );
#else
    testHierarchisationDehierarchisation(grid, level, &parabola, 1e-3 );
#endif
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationFundamentalSpline) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    for (int degree = 1; degree < 6; degree += 2) {
      Grid* grid = Grid::createFundamentalSplineGrid(dim, degree);
#if USE_DOUBLE_PRECISION
      testHierarchisationDehierarchisation(grid, level, &parabola, 1e-9, true);
#else
      testHierarchisationDehierarchisation(grid, level, &parabola, 0.01, true);
#endif
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationModFundamentalSpline) {
  int level = 5;

  for (int dim = 1; dim < 4; dim++) {
    for (int degree = 1; degree < 6; degree += 2) {
      Grid* grid = Grid::createModFundamentalSplineGrid(dim, degree);
#if USE_DOUBLE_PRECISION
      testHierarchisationDehierarchisation(grid, level, &parabola, 1e-9, true);
#else
      testHierarchisationDehierarchisation(grid, level, &parabola, 0.01, true);
#endif
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPoly) {
  int level = 5;
  int degree[4] = { 2, 3, 5, 8 };

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      Grid* grid = Grid::createPolyGrid(dim, degree[i]);
      testHierarchisationDehierarchisation(grid, level, &parabola, 0.0, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationPolyTruncatedBoundary) {
  int level = 5;
  int degree[4] = { 2, 3, 5, 8 };

  for (int dim = 1; dim < 4; dim++) {
    for (int i = 0; i < 4; i++) {
      Grid* grid = Grid::createPolyTruncatedBoundaryGrid(dim, degree[i]);
      testHierarchisationDehierarchisation(grid, level, &parabola, 0.0, true);
    }
  }
}

BOOST_AUTO_TEST_CASE(testHierarchisationStretchedTruncatedBoundary1D) {
  int dim = 1;
  int level = 5;
  Stretching1D* str1d = new Stretching1D();
  str1d->type = "log";
  str1d->x_0 = 0.;
  str1d->xsi = 10;
  DimensionBoundary* dimBound = new DimensionBoundary();
  dimBound->leftBoundary = 0.00001;
  dimBound->rightBoundary = 1.;
  Stretching stretch(dim, dimBound, str1d);
  Grid* grid = Grid::createLinearStretchedTruncatedBoundaryGrid(dim);
  grid->getStorage()->setStretching(stretch);
#if USE_DOUBLE_PRECISION == 1
  testHierarchisationDehierarchisation(grid, level, &parabolaBoundary, 1e-13, false, true);
#else
  testHierarchisationDehierarchisation(grid, level, &parabolaBoundary, 1e-4, false, true);
#endif
  delete str1d;
  delete dimBound;
}

BOOST_AUTO_TEST_CASE(testHierarchisationStretchedTruncatedBoundary3D) {
  int dim = 3;
  int level = 5;
  Stretching1D str1d;// = new Stretching1D();
  str1d.type = "sinh";
  str1d.x_0 = 1.;
  str1d.xsi = 10;
  DimensionBoundary dimBound;// = new DimensionBoundary();
  dimBound.leftBoundary = 0.001;
  dimBound.rightBoundary = 1.;

  std::vector<Stretching1D> stretch_vec;
  stretch_vec.push_back(str1d);
  stretch_vec.push_back(str1d);
  stretch_vec.push_back(str1d);
  std::vector<DimensionBoundary> dimBound_vec;
  dimBound_vec.push_back(dimBound);
  dimBound_vec.push_back(dimBound);
  dimBound_vec.push_back(dimBound);

  Stretching stretch(dim, dimBound_vec, stretch_vec);
  Grid* grid = Grid::createLinearStretchedTruncatedBoundaryGrid(dim);
  grid->getStorage()->setStretching(stretch);
#if USE_DOUBLE_PRECISION == 1
  testHierarchisationDehierarchisation(grid, level, &parabolaBoundary, 1e-12, false, true);
#else
  testHierarchisationDehierarchisation(grid, level, &parabolaBoundary, 1e-4, false, true);
#endif
}

BOOST_AUTO_TEST_SUITE_END()
