// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "GridCreator.hpp"

#include <vector>

void createSupportedGrids(size_t d, size_t p,
                          std::vector<std::unique_ptr<sgpp::base::Grid>>& grids) {
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createBsplineBoundaryGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createModBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createModBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createModFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createLinearBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createLinearClenshawCurtisGrid(d))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createModLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createWaveletGrid(d))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createWaveletBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<sgpp::base::Grid>(
                              sgpp::base::Grid::createModWaveletGrid(d))));
}

void createSampleGrid(sgpp::base::Grid& grid, size_t l, sgpp::optimization::ScalarFunction& f,
                      sgpp::base::DataVector& functionValues) {
  sgpp::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  // generate regular sparse grid
  gridStorage.emptyStorage();
  grid.getGenerator().regular(l);
  const size_t n = gridStorage.getSize();
  sgpp::base::DataVector x(d);

  functionValues.resize(n);

  for (size_t i = 0; i < n; i++) {
    sgpp::base::GridIndex& gp = *gridStorage[i];

    // don't forget to set the point distribution to Clenshaw-Curtis
    // if necessary (currently not done automatically)
    if (grid.getType() == sgpp::base::GridType::BsplineClenshawCurtis ||
        grid.getType() == sgpp::base::GridType::ModBsplineClenshawCurtis ||
        grid.getType() == sgpp::base::GridType::LinearClenshawCurtis) {
      gp.setPointDistribution(
        sgpp::base::GridIndex::PointDistribution::ClenshawCurtis);
    }

    for (size_t t = 0; t < d; t++) {
      x[t] = gp.getCoord(t);
    }

    functionValues[i] = f.eval(x);
  }
}

void createSampleGrid(sgpp::base::Grid& grid, size_t l, sgpp::optimization::VectorFunction& f,
                      sgpp::base::DataMatrix& functionValues) {
  sgpp::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();
  const size_t m = f.getNumberOfComponents();

  // generate regular sparse grid
  gridStorage.emptyStorage();
  grid.getGenerator().regular(l);
  const size_t n = gridStorage.getSize();
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector fx(m);

  functionValues.resize(n, m);

  for (size_t i = 0; i < n; i++) {
    sgpp::base::GridIndex& gp = *gridStorage[i];

    // don't forget to set the point distribution to Clenshaw-Curtis
    // if necessary (currently not done automatically)
    if (grid.getType() == sgpp::base::GridType::BsplineClenshawCurtis ||
        grid.getType() == sgpp::base::GridType::ModBsplineClenshawCurtis ||
        grid.getType() == sgpp::base::GridType::LinearClenshawCurtis) {
      gp.setPointDistribution(
        sgpp::base::GridIndex::PointDistribution::ClenshawCurtis);
    }

    for (size_t t = 0; t < d; t++) {
      x[t] = gp.getCoord(t);
    }

    f.eval(x, fx);
    functionValues.setRow(i, fx);
  }
}
