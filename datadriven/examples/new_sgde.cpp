// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/application/SparseGridDensityEstimator.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp_optimization.hpp>

#include <functional>
#include <random>
#include <vector>
#ifdef USE_CGAL
#include <CGAL/MP_Float.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_solution.h>
#include <CGAL/basic.h>
typedef CGAL::MP_Float ET;
typedef CGAL::Quadratic_program_solution<ET> Solution;
typedef CGAL::Quadratic_program_from_iterators<
    double**,                                              // for A
    double*,                                               // for b
    CGAL::Const_oneset_iterator<CGAL::Comparison_result>,  // for r
    bool*,                                                 // for fl
    double*,                                               // for l
    bool*,                                                 // for fu
    double*,                                               // for u
    double**,                                              // for D
    double*>                                               // for c
    Program;
#endif

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::GridType;
using sgpp::base::OperationEval;
using sgpp::base::OperationMatrix;
using sgpp::base::OperationMultipleEval;
using sgpp::datadriven::DensitySystemMatrix;

void calc_residual(Grid& grid, DataMatrix& samples, double lambda, DataVector& alpha,
                   DataVector& result) {
  size_t numSamples = samples.getNrows();
  size_t storage_size = grid.getStorage().getSize();
  OperationMatrix* A_op = sgpp::op_factory::createOperationLTwoDotProduct(grid);
  OperationMultipleEval* B_op = sgpp::op_factory::createOperationMultipleEval(grid, samples);
  OperationMatrix* C_op = sgpp::op_factory::createOperationLaplace(grid);
  DensitySystemMatrix AlambC(A_op, B_op, C_op, lambda, numSamples);
  DataVector b(storage_size);
  AlambC.generateb(b);
  AlambC.mult(alpha, result);
  result.sub(b);
}

void randu(DataVector& rvar, std::mt19937& generator) {
  std::normal_distribution<double> distribution(0.5, 0.1);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
    // make sure the sample is in the unit cube
    while (rvar[j] < 0 || rvar[j] > 1) rvar[j] = distribution(generator);
  }
}

void createSamples(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();
  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

void checkPositive(Grid& grid, const DataVector& alpha) {
  GridStorage* gridStorage = &grid.getStorage();
  size_t storage_size = gridStorage->getSize();
  size_t dims = grid.getDimension();
  GridPoint gp;
  DataVector coords(dims);
  OperationEval* op_eval = sgpp::op_factory::createOperationEval(grid);
  for (size_t i = 0; i < storage_size; i++) {
    gp = gridStorage->getPoint(i);
    gp.getStandardCoordinates(coords);
    if (op_eval->eval(alpha, coords) < 0) {
      std::cout << "f(" << coords.toString() << ")=" << op_eval->eval(alpha, coords) << std::endl;
    }
  }
}

#ifdef USE_CGAL
/**
 * solve new sparse grid density estimation using CGAL
 * @param grid
 * @param samples sample data
 * @param lambda for calculating the system matrix
 * @param result the best alpha vector found by the optimizer
 */
void solve_cgal(Grid& grid, DataMatrix& samples, double lambda, DataVector& result) {
  // CGAL documentation: https://doc.cgal.org/latest/QP_solver/index.html
  GridStorage* gridStorage = &(grid.getStorage());
  size_t storage_size = gridStorage->getSize();
  std::cout << "storage size:" << storage_size << std::endl;
  size_t numSamples = samples.getNrows();
  size_t dims = samples.getNcols();
  DataMatrix M(storage_size, storage_size);
  DataMatrix C(storage_size, storage_size);
  // loading M matrix
  OperationMatrix* A_op = sgpp::op_factory::createOperationLTwoDotExplicit(&M, grid);
  OperationMultipleEval* B_op = sgpp::op_factory::createOperationMultipleEval(grid, samples);
  // loading C matrix
  OperationMatrix* C_op = sgpp::op_factory::createOperationLaplaceExplicit(&C, grid);
  DataVector b(storage_size);
  DataVector q(storage_size);
  DensitySystemMatrix AlambC(A_op, B_op, C_op, lambda, numSamples);
  AlambC.generateb(b);

  // setting up System matrix: M + lambda*C
  C.mult(lambda);
  M.add(C);

  // Skalarproduct term of quadratic program: q.T x
  // where q.T = M.T b.T
  M.mult(b, q);
  q.mult(-1);

  // Quadratic program matrix P = M*M.T (M is symmetric)
  double** P_it = new double*[storage_size];
  for (size_t i = 0; i < storage_size; i++) {
    P_it[i] = new double[storage_size];
  }
  // multiplication P = M*M.T
  for (size_t i = 0; i < storage_size; i++) {
    DataVector col(storage_size);
    DataVector tmp(storage_size);
    M.getColumn(i, col);
    M.mult(col, tmp);
    P_it[i] = new double[storage_size];
    for (size_t j = 0; j <= i; j++) {
      P_it[i][j] = tmp.get(j);
      P_it[j][i] = tmp.get(j);
    }
  }

  // ---printing
  for (size_t i = 0; i < storage_size; i++) {
    for (size_t j = 0; j < storage_size; j++) {
      std::cout << P_it[j][i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "---------------" << std::endl;
  // ---printing

  // getting all grid points for interpolation matrix
  DataMatrix gridPoints(storage_size, dims);
  GridPoint gp;
  for (size_t i = 0; i < storage_size; i++) {
    gp = gridStorage->getPoint(i);
    DataVector coords(dims);
    gp.getStandardCoordinates(coords);
    gridPoints.setRow(i, coords);
  }

  // interpolation matrix (grid point values when multiplied with alpha-vec)
  OperationMultipleEval* G_op = sgpp::op_factory::createOperationMultipleEval(grid, gridPoints);
  double** G_it = new double*[storage_size];
  for (size_t i = 0; i < storage_size; i++) {
    G_it[i] = new double[storage_size];
  }
  // G_it[i] should contain the i-th COLUMN of G
  for (size_t i = 0; i < storage_size; i++) {
    DataVector tmp_result(storage_size);
    DataVector alpha(storage_size, 0.0);
    alpha.set(i, 1.0);
    G_op->mult(alpha, tmp_result);
    for (size_t j = 0; j < storage_size; j++) {
      G_it[i][j] = tmp_result.get(j);
    }
  }

  // ---printing
  for (size_t i = 0; i < storage_size; i++) {
    for (size_t j = 0; j < storage_size; j++) {
      std::cout << G_it[j][i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << q.toString() << std::endl;
  // ---printing

  // constraint relation (i.e. greater than zero)
  CGAL::Const_oneset_iterator<CGAL::Comparison_result> r(CGAL::LARGER);

  // lower and upper bounds for alpha are unused:
  bool* bounded = new bool[storage_size];
  for (size_t i = 0; i < storage_size; i++) {
    bounded[i] = false;
  }
  DataVector bounds(storage_size, 0.0);
  // define the quadratic Programm
  Program qp(static_cast<int>(storage_size), static_cast<int>(storage_size),  // size of problem
             // constraints
             G_it, bounds.getPointer(), r,
             // bounds
             bounded, bounds.getPointer(), bounded, bounds.getPointer(),
             // optimization goal
             P_it, q.getPointer());
  Solution s = CGAL::solve_quadratic_program(qp, ET());
  Solution::Variable_value_iterator it = s.variable_values_begin();
  for (size_t i = 0; i < storage_size; i++) {
    result.set(i, to_double(*it));
    it++;
  }
  std::cout << "objective function:" << to_double(s.objective_value()) << std::endl;
  // free up
  for (size_t i = 0; i < storage_size; i++) {
    delete P_it[i];
    delete G_it[i];
  }
  delete P_it;
  delete G_it;
  delete bounded;
  // delete A_op;
  // delete B_op;
  // delete C_op;
  // delete G_op;
}
#endif

int main(int argc, char** argv) {
  size_t d = 2;
  int level = 3;
  GridType gridType = sgpp::base::GridType::Linear;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = d;
  gridConfig.level_ = level;
  gridConfig.type_ = gridType;
  gridConfig.maxDegree_ = 5;
  DataMatrix trainSamples(1000, d);
  createSamples(trainSamples, 1234567);
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createGrid(gridConfig));
  grid->getGenerator().regular(gridConfig.level_);
  size_t storage_size = grid->getStorage().getSize();
  DataVector result(storage_size);

#ifdef USE_CGAL
  double lambda = 0.0;
  solve_cgal(*grid, trainSamples, lambda, result);
  std::cout << "Best alpha:" << result.toString() << std::endl;
  DataVector residual(storage_size);
  calc_residual(*grid, trainSamples, lambda, result, residual);
  std::cout << "residual:" << residual.l2Norm() << std::endl;
  double al[17] = {18.90769281, -9.09911171, -9.0176581,  0.10785147,  -1.18044529, -0.28708884,
                   0.15533845,  -9.38380473, -9.00915268, -0.02802101, -2.24950203, -1.58041974,
                   -0.09580067, 5.12791209,  4.54942521,  4.66378412,  4.32550775};
  DataVector cmp_alpha(al, 17);
  // calc_residual(*grid, trainSamples, lambda, cmp_alpha, residual);
  std::cout << "compare residual:" << residual.l2Norm() << std::endl;
  // feasibility checking
  std::cout << "Non positive best_alpha:" << std::endl;
  // checkPositive(*grid, result);
  std::cout << "Non positive cmp_alpha:" << std::endl;
// checkPositive(*grid, cmp_alpha);
#endif
  return 0;
}
