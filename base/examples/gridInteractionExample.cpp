// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstdlib>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>

#include <vector>
#include <unordered_set>
#include <tuple>

void decodeCoords(sgpp::base::DataVector& coords, std::vector<bool>& result) {
    for (size_t i = 0; i < coords.getSize(); ++i) {
        result[i] = coords[i] != 0.5;
    }
}


int main(int argc, char **argv) {
  auto dimensions = 3;
  auto level = 5;

  auto terms = std::unordered_set<std::vector<bool>>();
  auto realCoords = std::unordered_set<std::vector<bool>>();
  auto coords = sgpp::base::DataVector(dimensions);
  auto boolCoords = std::vector<bool>(dimensions);
  bool isInserted = false;

  // Add bias.
  terms.insert(std::vector<bool>(dimensions, false));


  // Add all variables, without interaction
  for (auto dim = 0; dim < dimensions; ++dim) {
      auto vec = std::vector<bool>(dimensions, false);
      vec[dim] = true;
      terms.insert(vec);
  }


  {
  auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
  auto& storage = grid->getStorage();
  auto generator = sgpp::base::HashGenerator();

  generator.regular_inter(storage, level, terms);

  std::cout << "Grid size with one level terms is " << grid->getSize() << std::endl;
  for (size_t i = 0; i < grid->getSize(); ++i) {
    sgpp::base::GridStorage::index_pointer gridIndex = storage.get(i);
    gridIndex->getCoords(coords);
    decodeCoords(coords, boolCoords);
    std::tie(std::ignore, isInserted) = realCoords.insert(boolCoords);
    if (isInserted) {
        std::cout << "New point" << coords.toString() << std::endl;
    }
  }
  }

  //Add all two-level-interactions
  for (auto i = 0; i < dimensions; ++i) {
      for (auto j = 0; j < dimensions; ++j) {
          auto vec = std::vector<bool>(dimensions, false);
          vec[i] = true;
          vec[j] = true;
          terms.insert(vec);
      }
  }

  {
  auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
  auto& storage = grid->getStorage();
  auto generator = sgpp::base::HashGenerator();

  generator.regular_inter(storage, level, terms);

  std::cout << "Grid size with two level terms is " << grid->getSize() << std::endl;
  for (size_t i = 0; i < grid->getSize(); ++i) {
    sgpp::base::GridStorage::index_pointer gridIndex = storage.get(i);
    gridIndex->getCoords(coords);
    decodeCoords(coords, boolCoords);
    std::tie(std::ignore, isInserted) = realCoords.insert(boolCoords);
    if (isInserted) {
        std::cout << "New point" << coords.toString() << std::endl;
    }
  }
  }

  //Add three-level-interaction
  auto vec = std::vector<bool>(dimensions, true);
  terms.insert(vec);

  {
  auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
  auto& storage = grid->getStorage();
  auto generator = sgpp::base::HashGenerator();

  generator.regular_inter(storage, level, terms);

  std::cout << "Grid size with three level terms is " << grid->getSize() << std::endl;

  for (size_t i = 0; i < grid->getSize(); ++i) {
    sgpp::base::GridStorage::index_pointer gridIndex = storage.get(i);
    gridIndex->getCoords(coords);
    decodeCoords(coords, boolCoords);
    std::tie(std::ignore, isInserted) = realCoords.insert(boolCoords);
    if (isInserted) {
        std::cout << "New point" << coords.toString() << std::endl;
    }
  }
  }
}



