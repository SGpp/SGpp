// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "GridCreator.hpp"

#include <vector>

void createSupportedGrids(size_t d, size_t p,
                          std::vector<std::unique_ptr<SGPP::base::Grid>>& grids) {
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createBsplineBoundaryGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createModBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createModBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createModFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createLinearBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createLinearClenshawCurtisGrid(d))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createModLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createWaveletGrid(d))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createWaveletBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<SGPP::base::Grid>(
                              SGPP::base::Grid::createModWaveletGrid(d))));
}

void createSampleGrid(SGPP::base::Grid& grid, size_t l, SGPP::optimization::ScalarFunction& f,
                      SGPP::base::DataVector& functionValues) {
  SGPP::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  // generate regular sparse grid
  gridStorage.emptyStorage();
  std::unique_ptr<SGPP::base::GridGenerator>(grid.createGridGenerator())->regular(l);
  const size_t n = gridStorage.size();
  SGPP::base::DataVector x(d);

  functionValues.resize(n);

  for (size_t i = 0; i < n; i++) {
    SGPP::base::GridIndex& gp = *gridStorage[i];

    // don't forget to set the point distribution to Clenshaw-Curtis
    // if necessary (currently not done automatically)
    if (grid.getType() == SGPP::base::GridType::BsplineClenshawCurtis ||
        grid.getType() == SGPP::base::GridType::ModBsplineClenshawCurtis ||
        grid.getType() == SGPP::base::GridType::LinearClenshawCurtis) {
      gp.setPointDistribution(
        SGPP::base::GridIndex::PointDistribution::ClenshawCurtis);
    }

    for (size_t t = 0; t < d; t++) {
      x[t] = gp.getCoord(t);
    }

    functionValues[i] = f.eval(x);
  }
}

void createSampleGrid(SGPP::base::Grid& grid, size_t l, SGPP::optimization::VectorFunction& f,
                      SGPP::base::DataMatrix& functionValues) {
  SGPP::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();
  const size_t m = f.getNumberOfComponents();

  // generate regular sparse grid
  gridStorage.emptyStorage();
  std::unique_ptr<SGPP::base::GridGenerator>(grid.createGridGenerator())->regular(l);
  const size_t n = gridStorage.size();
  SGPP::base::DataVector x(d);
  SGPP::base::DataVector fx(m);

  functionValues.resize(n, m);

  for (size_t i = 0; i < n; i++) {
    SGPP::base::GridIndex& gp = *gridStorage[i];

    // don't forget to set the point distribution to Clenshaw-Curtis
    // if necessary (currently not done automatically)
    if (grid.getType() == SGPP::base::GridType::BsplineClenshawCurtis ||
        grid.getType() == SGPP::base::GridType::ModBsplineClenshawCurtis ||
        grid.getType() == SGPP::base::GridType::LinearClenshawCurtis) {
      gp.setPointDistribution(
        SGPP::base::GridIndex::PointDistribution::ClenshawCurtis);
    }

    for (size_t t = 0; t < d; t++) {
      x[t] = gp.getCoord(t);
    }

    f.eval(x, fx);
    functionValues.setRow(i, fx);
  }
}
