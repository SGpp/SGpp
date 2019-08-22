// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::OperationMultipleEval;

BOOST_AUTO_TEST_SUITE(TestOperationMultipleEvalInterModLinear)

// zero alpha case -> evaluation must return alpha
BOOST_AUTO_TEST_CASE(zeroAlpha) {
  srand(static_cast<unsigned>(time(nullptr)));
  const size_t dim = size_t(4+10.*rand()/static_cast<double>(RAND_MAX));
  std::vector<std::vector<size_t>> interactions = std::vector<std::vector<size_t>>();

  std::vector<size_t> tmp = std::vector<size_t> ();

  // add empty interaction
  interactions.push_back(tmp);

  // add all unit interactions
  for (size_t i = 0; i < dim; i++) {
    tmp = std::vector<size_t> ();
    tmp.push_back(i);
    interactions.push_back(tmp);
  }

  // add ~50% of all binary interactions
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = i+1; j < dim; j++) {
      if (rand()/static_cast<double>(RAND_MAX) > .5) {
        tmp = std::vector<size_t> ();
        tmp.push_back(i);
        tmp.push_back(j);
        interactions.push_back(tmp);
      }
    }
  }

  std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
  grid->getGenerator().regularInter(3, interactions, 0.);

  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();

  DataVector alpha(N);

  // set all alphas to zero
  for (int i = 0; i < static_cast<int>(N); ++i) {
    alpha[i] = static_cast<double>(0.);
  }

  // generate 20 random points
  const size_t numberDataPoints = 20;
  double **points = new double* [numberDataPoints];

  for (size_t i = 0; i < numberDataPoints; i++) {
    points[i] = new double[dim];
    for (size_t j = 0; j < dim; j++) {
      points[i][j] = rand()/static_cast<double>(RAND_MAX);
    }
  }

  DataVector result(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result[i] = static_cast<double>(0);
  }

  DataMatrix dataset(numberDataPoints, dim);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    DataVector temp(dim);

    for (int j = 0; j < static_cast<int>(dim); ++j) {
      temp[j] = points[i][j];
    }

    dataset.setRow(i, temp);
  }

  sgpp::op_factory::createOperationMultipleEvalInter(*grid, dataset, interactions)
  ->mult(alpha, result);

  BOOST_TEST_MESSAGE(alpha.toString() + "\n");
  BOOST_TEST_MESSAGE(result.toString() + "\n");

  for (size_t i = 0; i < numberDataPoints; i++) {
    BOOST_CHECK_CLOSE(result[i], 0, 1e-7);
  }

  for (size_t i = 0; i < numberDataPoints; i++) {
    delete[] points[i];
  }
  delete[] points;
}

// when all coefficients have the same value the grid
// should be symmetric in all dimensions
BOOST_AUTO_TEST_CASE(symmetry) {
  srand(static_cast<unsigned>(time(nullptr)));
  const size_t dim = size_t(4+10*(rand()/static_cast<double>(RAND_MAX)));
  std::vector<std::vector<size_t>> interactions = std::vector<std::vector<size_t>>();

  std::vector<size_t> tmp = std::vector<size_t> ();

  // add empty interaction
  interactions.push_back(tmp);

  // add all unit interactions
  for (size_t i = 0; i < dim; i++) {
    tmp = std::vector<size_t> ();
    tmp.push_back(i);
    interactions.push_back(tmp);
  }

  // add ~50% of all binary interactions
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = i+1; j < dim; j++) {
      if (rand()/static_cast<double>(RAND_MAX) > .5) {
        tmp = std::vector<size_t> ();
        tmp.push_back(i);
        tmp.push_back(j);
        interactions.push_back(tmp);
      }
    }
  }

  std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
  grid->getGenerator().regularInter(3, interactions, 0.);

  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();

  DataVector alpha(N);

  double coef = 1+rand()/static_cast<double>(RAND_MAX);
  // set all alphas to the same random value
  for (int i = 0; i < static_cast<int>(N); ++i) {
    alpha[i] = static_cast<double>(coef);
  }

  // generate 20 random points
  const size_t numberDataPoints = 20;
  double **points = new double* [numberDataPoints];

  for (size_t i = 0; i < numberDataPoints; i++) {
    points[i] = new double[dim];
    for (size_t j = 0; j < dim; j++) {
      points[i][j] = rand()/static_cast<double>(RAND_MAX);
    }
  }

  DataVector result(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result[i] = static_cast<double>(0);
  }

  DataMatrix dataset(numberDataPoints, dim);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    DataVector temp(dim);

    for (int j = 0; j < static_cast<int>(dim); ++j) {
      temp[j] = points[i][j];
    }

    dataset.setRow(i, temp);
  }

  DataVector result_mirrored(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result_mirrored[i] = static_cast<double>(0);
  }
  // same dataset with randomly mirrored dimensions
  DataMatrix dataset_mirrored(numberDataPoints, dim);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    DataVector temp(dim);

    for (int j = 0; j < static_cast<int>(dim); ++j) {
      if (rand()/static_cast<double>(RAND_MAX) < 0.25)
        temp[j] = 1. - points[i][j];
      else
        temp[j] = points[i][j];
    }

    dataset_mirrored.setRow(i, temp);
  }

  sgpp::op_factory::createOperationMultipleEvalInter(*grid, dataset, interactions)
  ->mult(alpha, result);


  sgpp::op_factory::createOperationMultipleEvalInter(*grid, dataset_mirrored, interactions)
  ->mult(alpha, result_mirrored);

  BOOST_TEST_MESSAGE(result_mirrored.toString() + "\n");
  BOOST_TEST_MESSAGE(result.toString() + "\n");

  for (size_t i = 0; i < numberDataPoints; i++) {
    BOOST_CHECK_CLOSE(result[i], result_mirrored[i], 1e-7);
  }

  for (size_t i = 0; i < numberDataPoints; i++) {
    delete[] points[i];
  }
  delete[] points;
}

// check if the interaction evaluation returns
// the same result as the regular one
BOOST_AUTO_TEST_CASE(regularEvaluation) {
  srand(static_cast<unsigned>(time(nullptr)));
  const size_t dim = size_t(4+10*(rand()/static_cast<double>(RAND_MAX)));
  std::vector<std::vector<size_t>> interactions = std::vector<std::vector<size_t>>();

  std::vector<size_t> tmp = std::vector<size_t> ();

  // add empty interaction
  interactions.push_back(tmp);

  // add all unit interactions
  for (size_t i = 0; i < dim; i++) {
    tmp = std::vector<size_t> ();
    tmp.push_back(i);
    interactions.push_back(tmp);
  }

  // add ~50% of all binary interactions
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = i+1; j < dim; j++) {
      if (rand()/static_cast<double>(RAND_MAX) > .5) {
        tmp = std::vector<size_t> ();
        tmp.push_back(i);
        tmp.push_back(j);
        interactions.push_back(tmp);
      }
    }
  }

  std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
  grid->getGenerator().regularInter(4, interactions, 0.);

  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();

  DataVector alpha(N);


  double coef = 1+rand()/static_cast<double>(RAND_MAX);
  // set all alphas to the same random value
  for (int i = 0; i < static_cast<int>(N); ++i) {
    alpha[i] = static_cast<double>(coef);
  }

  // generate 20 random points
  const size_t numberDataPoints = 20;
  double **points = new double* [numberDataPoints];

  for (size_t i = 0; i < numberDataPoints; i++) {
    points[i] = new double[dim];
    for (size_t j = 0; j < dim; j++) {
      points[i][j] = rand()/static_cast<double>(RAND_MAX);
    }
  }

  DataVector result(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result[i] = static_cast<double>(0);
  }

  DataMatrix dataset(numberDataPoints, dim);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    DataVector temp(dim);

    for (int j = 0; j < static_cast<int>(dim); ++j) {
      temp[j] = points[i][j];
    }

    dataset.setRow(i, temp);
  }

  DataVector result_reg(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result_reg[i] = static_cast<double>(0);
  }

  sgpp::op_factory::createOperationMultipleEvalInter(*grid, dataset, interactions)
  ->mult(alpha, result);


  sgpp::op_factory::createOperationMultipleEval(*grid, dataset)
  ->mult(alpha, result_reg);

  BOOST_TEST_MESSAGE(result_reg.toString() + "\n");
  BOOST_TEST_MESSAGE(result.toString() + "\n");

  for (size_t i = 0; i < numberDataPoints; i++) {
    BOOST_CHECK_CLOSE(result[i], result_reg[i], 1e-7);
  }

  for (size_t i = 0; i < numberDataPoints; i++) {
    delete[] points[i];
  }
  delete[] points;
}

BOOST_AUTO_TEST_SUITE_END()
