// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "GridCreator.hpp"

#include <vector>


void createSupportedGrids(
    size_t d, size_t p, std::vector<std::unique_ptr<sgpp::base::Grid>>& grids) {
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createBsplineClenshawCurtisGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModBsplineClenshawCurtisGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createFundamentalNakSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createFundamentalSplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createFundamentalSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModFundamentalSplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createWeaklyFundamentalNakSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createWeaklyFundamentalSplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModWeaklyFundamentalNakSplineGrid(d, p)));
  grids.push_back(
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createLinearBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createLinearClenshawCurtisBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createLinearClenshawCurtisGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModLinearGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModNakBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createNaturalBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createNakBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createWaveletGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createWaveletBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModWaveletGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
     sgpp::base::Grid::createNakBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createNakBsplineExtendedGrid(d, p)));
}

void createSampleGrid(sgpp::base::Grid& grid, size_t l,
                      sgpp::base::ScalarFunction& f,
                      sgpp::base::DataVector& functionValues) {
  sgpp::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  // generate regular sparse grid
  gridStorage.clear();
  grid.getGenerator().regular(l);
  const size_t n = gridStorage.getSize();
  sgpp::base::DataVector x(d);

  functionValues.resize(n);

  for (size_t i = 0; i < n; i++) {
    x = gridStorage.getCoordinates(gridStorage[i]);
    functionValues[i] = f.eval(x);
  }
}

void createSampleGrid(sgpp::base::Grid& grid, size_t l,
                      sgpp::base::VectorFunction& f,
                      sgpp::base::DataMatrix& functionValues) {
  sgpp::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();
  const size_t m = f.getNumberOfComponents();

  // generate regular sparse grid
  gridStorage.clear();
  grid.getGenerator().regular(l);
  const size_t n = gridStorage.getSize();
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector fx(m);

  functionValues.resize(n, m);

  for (size_t i = 0; i < n; i++) {
    x = gridStorage.getCoordinates(gridStorage[i]);
    f.eval(x, fx);
    functionValues.setRow(i, fx);
  }
}
