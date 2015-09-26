#include "GridCreator.hpp"

using namespace SGPP;

void createSupportedGrids(size_t d, size_t p,
                          std::vector<std::unique_ptr<base::Grid>>& grids) {
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createBsplineBoundaryGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createLinearBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createLinearClenshawCurtisGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createWaveletGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createWaveletBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModWaveletGrid(d))));
}

void createSampleGrid(base::Grid& grid, size_t l, ScalarFunction& f,
                      base::DataVector& functionValues) {
  base::GridStorage& gridStorage = *grid.getStorage();
  const size_t d = gridStorage.dim();

  // generate regular sparse grid
  gridStorage.emptyStorage();
  std::unique_ptr<base::GridGenerator> gridGen(grid.createGridGenerator());
  gridGen->regular(l);
  const size_t n = gridStorage.size();
  base::DataVector x(d);

  functionValues.resize(n);

  for (size_t i = 0; i < n; i++) {
    base::GridIndex& gp = *gridStorage[i];

    // don't forget to set the point distribution to Clenshaw-Curtis
    // if necessary (currently not done automatically)
    if (grid.getType() == base::GridType::BsplineClenshawCurtis ||
        grid.getType() == base::GridType::ModBsplineClenshawCurtis ||
        grid.getType() == base::GridType::LinearClenshawCurtis) {
      gp.setPointDistribution(
        base::GridIndex::PointDistribution::ClenshawCurtis);
    }

    for (size_t t = 0; t < d; t++) {
      x[t] = gp.getCoord(t);
    }

    functionValues[i] = f.eval(x);
  }
}
