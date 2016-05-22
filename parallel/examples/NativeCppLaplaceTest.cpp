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
    std::cout << std::endl
              << "Usage " << argv[0]
              << " [dim] [level] [bound_left] [bound_right] [dirichlet] [M/V] [I/B] [Lambda(0/1)]"
              << std::endl
              << std::endl;
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

  sgpp::base::BoundingBox1D* myBoundaries = new sgpp::base::BoundingBox1D[dim];

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

  sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(dim, myBoundaries);
  delete[] myBoundaries;

  // Generate lambda for test if needed
  sgpp::base::DataVector lambda(dim);

  if (lmb == 1) {
    for (size_t d = 0; d < dim; d++)
      lambda.set(d, (static_cast<double>(d + 1) / static_cast<double>(dim)));
  } else {
    lambda.setAll(1.0);
  }

  // Test inner grid
  sgpp::base::Grid* myGrid;

  if (grid_selection == "I") {
    myGrid = new sgpp::base::LinearGrid(*myBoundingBox);
  } else if (grid_selection == "B") {
    myGrid = new sgpp::base::LinearBoundaryGrid(*myBoundingBox);
  } else {
    std::cout << std::endl
              << "Usage " << argv[0]
              << " [dim] [level] [bound_left] [bound_right] [dirichlet] [M/V] [I/B] [Lambda(0/1)]"
              << std::endl
              << std::endl;
    exit(-1);
  }

  myGrid->getGenerator().regular(level);

  sgpp::base::GridStorage* myGridStorage = &myGrid->getStorage();
  size_t gridsize = myGridStorage->getSize();
  sgpp::base::DataVector* alpha = new sgpp::base::DataVector(gridsize);
  sgpp::base::DataVector* result_updown = new sgpp::base::DataVector(gridsize);
  sgpp::base::DataVector* result_vector = new sgpp::base::DataVector(gridsize);

  std::unique_ptr<sgpp::base::OperationMatrix> updown;
  std::unique_ptr<sgpp::base::OperationMatrix> vect;

  if (lmb == 1) {
    updown = sgpp::op_factory::createOperationLaplace(*myGrid, lambda);
    vect = sgpp::op_factory::createOperationLaplaceVectorized(
        *myGrid, lambda,
        sgpp::parallel::X86SIMD);  /// @todo: check for parallelization type
  } else {
    updown = sgpp::op_factory::createOperationLaplace(*myGrid);
    vect = sgpp::op_factory::createOperationLaplaceVectorized(
        *myGrid,
        sgpp::parallel::X86SIMD);  /// @todo: check for parallelization type
  }

  std::cout << std::endl;
  std::cout << "Laplace Test Application" << std::endl;
  std::cout << "========================" << std::endl << std::endl;
  std::cout << "Gridtype:     " << static_cast<int>(myGrid->getType())
            << std::endl;  /// @todo:output for gridtype
  std::cout << "Dimensions:   " << dim << std::endl;
  std::cout << "Levels:       " << level << std::endl;
  std::cout << "Left Bound.:  " << bound_left << std::endl;
  std::cout << "Right Bound.: " << bound_right << std::endl;
  std::cout << "Gridsize:     " << gridsize << std::endl;
  std::cout << "chosen Test:  " << test << std::endl;
  std::cout << "Lambda:       ";

  for (size_t d = 0; d < dim; d++) std::cout << lambda.get(d) << " ";

  std::cout << std::endl << std::endl;

  double max_error = 0.0;

  if (test == "M") {
    alpha->setAll(0.0);
    sgpp::base::DataMatrix* matrix_updown = new sgpp::base::DataMatrix(gridsize, gridsize);
    sgpp::base::DataMatrix* matrix_vector = new sgpp::base::DataMatrix(gridsize, gridsize);

    for (size_t i = 0; i < gridsize; i++) {
      alpha->set(i, 1.0);
      result_updown->setAll(0.0);
      result_vector->setAll(0.0);

      updown->mult(*alpha, *result_updown);
      matrix_updown->setColumn(i, *result_updown);
      vect->mult(*alpha, *result_vector);
      matrix_vector->setColumn(i, *result_vector);

      // if ((i % (gridsize/20)) == 0) std::cout << "." ;

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
    std::cout << std::endl
              << "Usage " << argv[0]
              << " [dim] [level] [bound_left] [bound_right] [dirichlet] [M/V] [I/B]" << std::endl
              << std::endl;
    exit(-1);
  }

  std::cout << "Max error for vectorized version is " << max_error << std::endl << std::endl;

  delete result_vector;
  delete result_updown;
  delete alpha;
  delete myGrid;

  return 0;
}
