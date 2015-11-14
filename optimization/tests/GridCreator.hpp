#ifndef GRID_CREATOR_HPP
#define GRID_CREATOR_HPP

#include <vector>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>

using namespace SGPP;
using namespace SGPP::optimization;

void createSupportedGrids(size_t d, size_t p,
                          std::vector<std::unique_ptr<base::Grid>>& grids);

void createSampleGrid(base::Grid& grid, size_t l, ScalarFunction& f,
                      base::DataVector& functionValues);

void createSampleGrid(base::Grid& grid, size_t l, VectorFunction& f,
                      base::DataMatrix& functionValues);

#endif /* GRID_CREATOR_HPP */
