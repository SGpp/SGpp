// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRID_CREATOR_HPP
#define GRID_CREATOR_HPP

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <vector>

void createSupportedGrids(size_t d, size_t p,
                          std::vector<std::unique_ptr<sgpp::base::Grid>>& grids);

void createSampleGrid(sgpp::base::Grid& grid, size_t l, sgpp::base::ScalarFunction& f,
                      sgpp::base::DataVector& functionValues);

void createSampleGrid(sgpp::base::Grid& grid, size_t l, sgpp::base::VectorFunction& f,
                      sgpp::base::DataMatrix& functionValues);

#endif /* GRID_CREATOR_HPP */
