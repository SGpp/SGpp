// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <bits/stdc++.h>
#include <chrono>
#include <set>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::OperationMultipleEval;

void genAllInteractions(std::set<std::set<size_t>> &inter, size_t d) {
  std::set<size_t> tmp;
  inter.emplace(tmp);

  size_t count = static_cast<size_t>(pow(2., static_cast<double>(d)));
  for (size_t i = 1; i < count; ++i) {
    tmp.clear();
    for (size_t j = 0; j < d; ++j) {
      if ((i & (1 << j)) > 0) {
        tmp.emplace(j);
      }
    }
    inter.emplace(tmp);
  }
}

int main() {
  std::vector<size_t> dims = {1, 2, 3};
  size_t lvl_start = 3;
  size_t lvl_end = 5;

  for (auto const &d : dims) {
    std::set<std::set<size_t>> interactions;
    genAllInteractions(interactions, d);

    /*std::cout << std::endl << "~~~" << std::endl << std::endl;
    for (auto it1 = interactions.begin(); it1 != interactions.end(); ++it1) {
      std::cout << "[ ";
      for (auto it2 = it1->begin(); it2 != it1->end(); ++it2) {
        std::cout << *it2 << " ";
      }
      std::cout << "]" << std::endl;
    }
    std::cout << std::endl;*/

    for (size_t lvl = lvl_start; lvl <= lvl_end; ++lvl) {
      // get info on grid
      std::unique_ptr<Grid> grid(Grid::createModLinearGrid(d));
      grid->getGenerator().regularInter(lvl, interactions, 0.);
      GridStorage &gS = grid->getStorage();
      size_t N = gS.getSize();

      // Create alpha vector
      DataVector alpha(N);

      // set all alphas to zero
      for (int i = 0; i < static_cast<int>(N); ++i) {
        alpha[i] = static_cast<double>(0.);
      }

      // generate random points
      const size_t numberDataPoints = lvl * d * 10;
      double **points = new double *[numberDataPoints];

      for (size_t i = 0; i < numberDataPoints; i++) {
        points[i] = new double[d];
        for (size_t j = 0; j < d; j++) {
          points[i][j] = rand() / static_cast<double>(RAND_MAX);
        }
      }

      DataVector result(numberDataPoints);

      for (unsigned int i = 0; i < (numberDataPoints); ++i) {
        result[i] = static_cast<double>(0);
      }

      DataMatrix dataset(numberDataPoints, d);

      for (unsigned int i = 0; i < (numberDataPoints); ++i) {
        DataVector temp(d);

        for (int j = 0; j < static_cast<int>(d); ++j) {
          temp[j] = points[i][j];
        }

        dataset.setRow(i, temp);
      }

      // Make the runs; We use 10000 runs of each, to avoid measurement bias
      std::cout << "dim: " << d << ", lvl: " << lvl
                << ", no.of grid points: " << N
                << ", no. of data points: " << numberDataPoints << std::endl;

      auto start = std::chrono::high_resolution_clock::now();
      for (size_t k = 0; k < 10000; ++k) {
        // Run and measure time for interaction-based eval
        sgpp::op_factory::createOperationMultipleEvalInter(*grid, dataset,
                                                           interactions)
            ->mult(alpha, result);
      }
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      std::cout << "\tInteract 10000 runs\t: " << elapsed.count() << std::endl;

      // Create new grid
      grid = std::unique_ptr<Grid>(Grid::createModLinearGrid(d));
      grid->getGenerator().regular(lvl, 0.);

      start = std::chrono::high_resolution_clock::now();
      for (size_t k = 0; k < 10000; ++k) {
        // Run and measure time for regular eval
        sgpp::op_factory::createOperationMultipleEval(*grid, dataset)
            ->mult(alpha, result);
      }
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      std::cout << "\tStandard 10000 runs\t: " << elapsed.count() << std::endl;

      // Release memory
      for (size_t i = 0; i < numberDataPoints; i++) {
        delete[] points[i];
      }
      delete[] points;
    }
  }

  return 0;
}
