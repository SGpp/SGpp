#ifndef GRID_CREATOR_HPP
#define GRID_CREATOR_HPP

#include <vector>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/ObjectiveFunction.hpp>

using namespace SGPP;
using namespace SGPP::optimization;

void createSupportedGrids(size_t d, size_t p,
                          std::vector<std::unique_ptr<base::Grid>>& grids);

void createSampleGrid(base::Grid& grid, size_t l, ObjectiveFunction& f,
                      base::DataVector& functionValues);

#endif /* GRID_CREATOR_HPP */
