// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp_parallel.hpp>

#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>

/**
 * Testapplication for the Intel VTune Profiling Tool
 * and a measurement app for Sparse Grid Algorithms building blocks
 */
int main(int argc, char* argv[]) {
  size_t dim;
  int level;
  double bound_left;
  double bound_right;
  int dirichlet;
  int lmb;
  std::string test;
  std::string grid_selection;

  if (argc != 9) {
    std::cout << std::endl << "Usage " << argv[0] << " [dim] [level] [bound_left] [bound_right] [dirichlet] [M/V] [I/B] [Lambda(0/1)]" << std::endl << std::endl;
    exit(-1);
  }

  dim = atoi(argv[1]);
  level = atoi(argv[2]);
  bound_left = atof(argv[3]);
  bound_right = atof(argv[4]);
  dirichlet = atoi(argv[5]);
  test.assign(argv[6]);
  grid_selection.assign(argv[7]);
  lmb = atoi(argv[8]);

  SGPP::base::DimensionBoundary* myBoundaries = new SGPP::base::DimensionBoundary[dim];

  // set the bounding box
  for (size_t i = 0; i < dim; i++) {
    myBoundaries[i].leftBoundary = bound_left;
    myBoundaries[i].rightBoundary = bound_right;

    if (dirichlet == 0) {
      myBoundaries[i].bDirichletLeft = false;
      myBoundaries[i].bDirichletRight = false;
    } else {
      myBoundaries[i].bDirichletLeft = true;
      myBoundaries[i].bDirichletRight = true;
    }
  }

  SGPP::base::BoundingBox* myBoundingBox = new SGPP::base::BoundingBox(dim, myBoundaries);
  delete[] myBoundaries;

  // Generate lambda for test if needed
  SGPP::base::DataVector lambda(dim);

  if (lmb == 1) {
    for (size_t d = 0; d < dim; d++)
      lambda.set(d, ((double)(d + 1) / (double)dim));
  } else {
    lambda.setAll(1.0);
  }

  // Test inner grid
  SGPP::base::Grid* myGrid;

  if (grid_selection == "I") {
    myGrid = new SGPP::base::LinearGrid(*myBoundingBox);
  } else if (grid_selection == "B") {
    myGrid = new SGPP::base::LinearBoundaryGrid(*myBoundingBox);
  } else {
    std::cout << std::endl << "Usage " << argv[0] << " [dim] [level] [bound_left] [bound_right] [dirichlet] [M/V] [I/B] [Lambda(0/1)]" << std::endl << std::endl;
    exit(-1);
  }

  SGPP::base::GridGenerator* myGenerator = myGrid->createGridGenerator();
  myGenerator->regular(level);
  delete myGenerator;

  SGPP::base::GridStorage* myGridStorage = myGrid->getStorage();
  size_t gridsize = myGridStorage->size();
  SGPP::base::DataVector* alpha = new SGPP::base::DataVector(gridsize);
  SGPP::base::DataVector* result_updown = new SGPP::base::DataVector(gridsize);
  SGPP::base::DataVector* result_vector = new SGPP::base::DataVector(gridsize);

  SGPP::base::OperationMatrix* updown;
  SGPP::base::OperationMatrix* vect;

  if (lmb == 1) {
    updown = SGPP::op_factory::createOperationLaplace(*myGrid, lambda);
    vect = SGPP::op_factory::createOperationLaplaceVectorized(*myGrid, lambda);
  } else {
    updown = SGPP::op_factory::createOperationLaplace(*myGrid);
    vect = SGPP::op_factory::createOperationLaplaceVectorized(*myGrid);
  }

  std::cout << std::endl;
  std::cout << "Laplace Test Application" << std::endl;
  std::cout << "========================" << std::endl << std::endl;
  std::cout << "Gridtype:     " << myGrid->getType() << std::endl;
  std::cout << "Dimensions:   " << dim << std::endl;
  std::cout << "Levels:       " << level << std::endl;
  std::cout << "Left Bound.:  " << bound_left << std::endl;
  std::cout << "Right Bound.: " << bound_right << std::endl;
  std::cout << "Gridsize:     " << gridsize << std::endl;
  std::cout << "chosen Test:  " << test << std::endl;
  std::cout << "Lambda:       ";

  for (size_t d = 0; d < dim; d++)
    std::cout << lambda.get(d) << " ";

  std::cout << std::endl << std::endl;

  double max_error = 0.0;

  if (test == "M") {
    alpha->setAll(0.0);
    SGPP::base::DataMatrix* matrix_updown = new SGPP::base::DataMatrix(gridsize, gridsize);
    SGPP::base::DataMatrix* matrix_vector = new SGPP::base::DataMatrix(gridsize, gridsize);

    for (size_t i = 0; i < gridsize; i++) {
      alpha->set(i, 1.0);
      result_updown->setAll(0.0);
      result_vector->setAll(0.0);

      updown->mult(*alpha, *result_updown);
      matrix_updown->setColumn(i, *result_updown);
      vect->mult(*alpha, *result_vector);
      matrix_vector->setColumn(i, *result_vector);

      //if ((i % (gridsize/20)) == 0) std::cout << "." ;

      alpha->set(i, 0.0);
    }

    std::cout << std::endl;

    // compare result
    for (size_t i = 0; i < gridsize; i++) {
      for (size_t j = 0; j < gridsize; j++) {
        if (fabs(matrix_updown->get(i, j) - matrix_vector->get(i, j)) > max_error) {
          max_error = fabs(matrix_updown->get(i, j) - matrix_vector->get(i, j));
        }
      }
    }

#if 0
    std::cout << std::endl << std::endl;

    for (size_t i = 0; i < gridsize; i++) {
      for (size_t j = 0; j < gridsize; j++) {
        std::cout << matrix_updown->get(i, j) << " ";
      }

      std::cout << std::endl;
    }

    std::cout << std::endl << std::endl;

    for (size_t i = 0; i < gridsize; i++) {
      for (size_t j = 0; j < gridsize; j++) {
        std::cout << matrix_vector->get(i, j) << " ";
      }

      std::cout << std::endl;
    }

    std::cout << std::endl << std::endl;

    std::cout << myGrid->serialize() << std::endl;

    std::cout << std::endl << std::endl;
#endif
    delete matrix_updown;
    delete matrix_vector;
  } else if (test == "V") {
    alpha->setAll(1.0);
    updown->mult(*alpha, *result_updown);
    vect->mult(*alpha, *result_vector);

    // compare result
    for (size_t i = 0; i < gridsize; i++) {
      if (fabs(result_updown->get(i) - result_vector->get(i)) > max_error) {
        max_error = fabs(result_updown->get(i) - result_vector->get(i));
      }

    }
  } else {
    std::cout << std::endl << "Usage " << argv[0] << " [dim] [level] [bound_left] [bound_right] [dirichlet] [M/V] [I/B]" << std::endl << std::endl;
    exit(-1);
  }

  std::cout << "Max error for vectorized version is " << max_error << std::endl << std::endl;

  delete vect;
  delete updown;
  delete result_vector;
  delete result_updown;
  delete alpha;
  delete myGrid;

  return 0;
}