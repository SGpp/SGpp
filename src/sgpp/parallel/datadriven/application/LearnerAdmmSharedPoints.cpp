/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "parallel/datadriven/application/LearnerADMMSharedPoints.hpp"
#include "datadriven/algorithm/DMSystemMatrix.hpp"
#include "sgpp_mpi.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "base/exception/solver_exception.hpp"
#include "base/exception/application_exception.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <mpi.h>

//#include "omp.h"

#include "datadriven/algorithm/DMSystemMatrixExplicit.hpp"


using namespace sg::parallel;
using namespace sg::datadriven;
using namespace sg::base;

int LearnerADMMSharedPoints::Factorial(int nValue) {
  if (nValue < 0)
    throw new sg::base::application_exception("Factorial of negative number cannot be calculated");

  if (nValue == 0)
    return 1;

  int result = nValue--;

  for (; nValue > 1; nValue--)
    result *= nValue;

  return result;
}

int LearnerADMMSharedPoints::EvaluateBinomialCoefficient(int nValue, int nValue2) {
  int result;

  if (nValue2 == 1)
    return nValue;

  result = (Factorial(nValue)) / (Factorial(nValue2) * Factorial((nValue - nValue2)));
  nValue2 = result;
  return nValue2;
}

int LearnerADMMSharedPoints::levels2functions(int level) {
  int result = 0;
  int dim =  grid_->getStorage()->dim();

  for (int l = 0; l < level; l++) {
    result += (1 << l) * EvaluateBinomialCoefficient(l + dim - 1, dim - 1);
  }

  return result;
}

void LearnerADMMSharedPoints::createFactorizedSystem(DataMatrix& trainDataset, double lambda) {
  num_shared_functions = levels2functions(num_shared_levels);

  if (rank == 0)
    std::cout << "Number of shared points: " << num_shared_functions << std::endl;

  //double rho = 1.0;
  int grid_points = this->grid_->getSize();

  if (grid_points < num_shared_functions) {
    std::stringstream message("Number of shared points ");
    message << num_shared_functions << " is greater than number of grid points " << grid_points << std::endl;
    throw base::application_exception(message.str().c_str());
  }

  grid_points -= num_shared_functions;

  train_size = trainDataset.getNrows();
  int rows = train_size;

  // ceiling
  int part_size = (grid_points - num_shared_functions) / num_tasks + (((grid_points - num_shared_functions) % num_tasks) ? 1 : 0);
  first_index = part_size * rank + num_shared_functions;
  last_index = (part_size * (rank + 1) + num_shared_functions < grid_points) ? part_size * (rank + 1) - 1 + num_shared_functions : grid_points - 1;
  subproblem_length = last_index - first_index + 1  + num_shared_functions;


  //std::cout << "rank " << rank << " subproblem size " << subproblem_length << std::endl;
  if (rank == 0)
    std::cout << "building matrix" << std::endl;

  //Assemble matrix
  sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
  myStopwatch->start();

  DMSystemMatrixExplicit dmsystem(trainDataset, lambda);

  if (B_ != NULL)
    delete B_;

  // correct
  double* B1 = dmsystem.compileMatrix(this->grid_->getStorage(), 0, num_shared_functions - 1);
  double* B2 = dmsystem.compileMatrix(this->grid_->getStorage(), first_index, last_index);
  sleep(rank);

  B_ = new double[subproblem_length * train_size];

  for (int i = 0; i < train_size; i++) {
    for (int j = 0; j < num_shared_functions; j++) {
      B_[i * subproblem_length + j] = B1[i * num_shared_functions + j];
    }
  }

  for (int i = 0; i < train_size; i++) {
    for (int j = 0; j < subproblem_length - num_shared_functions ; j++) {
      B_[i * subproblem_length + j + num_shared_functions] = B2[i * (subproblem_length - num_shared_functions) + j ];
    }
  }


  if (system_matrix != NULL)
    delete[] system_matrix;

  system_matrix = new double[subproblem_length * subproblem_length];

  for (int i = 0; i < subproblem_length; i++) {
    for (int j = 0; j < subproblem_length; j++)
      system_matrix[i * subproblem_length + j] = (i == j) ? 1.0 : 0.0; //diagonal elements
  }

  // system_matrix = \rho*B^T*B + \lambda*system_matrix
  // TODO: should it be train_size*\lambda?
  //correct
  if (rank == 0)
    std::cout << "lambda value: " << lambda << std::endl;

  cblas_dsyrk(CblasRowMajor, CblasLower,
              CblasTrans, subproblem_length, rows, rho, B_,
              subproblem_length, lambda, system_matrix, subproblem_length);
  std::cout << "rank " << rank << " building matrix directly in " << myStopwatch->stop() << " seconds" << std::endl;

  //Cholesky factorization
  myStopwatch->start();

  int info;

  if (solution_method == Cholesky) {
    this->dpotrf("U", subproblem_length, system_matrix, subproblem_length, info);
    std::cout << "rank " << rank << " performing cholesky factorization in " << myStopwatch->stop() << " seconds" << std::endl;
  }

  delete myStopwatch;
  delete[] B1, B2;
}

void LearnerADMMSharedPoints::set_num_shared_levels(int num) {
  this->num_shared_levels = num;
}

int LearnerADMMSharedPoints::get_num_shared_levels() {
  return this->num_shared_levels;
}


LearnerTiming LearnerADMMSharedPoints::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
    const double lambda, double convergence_threshold) {
  LearnerTiming result;
  sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
  myStopwatch->start();

  // initializations
  this->createFactorizedSystem(trainDataset, lambda);
  int grid_size = this->grid_->getSize();
  int processes = num_tasks;
  double* Bixi = new double[train_size]; //(double*)std::calloc(train_size, sizeof(double));
  double* Bixi_shared = (double*)calloc(train_size, sizeof(double));
  double* Bixi_shared_sum = (double*)calloc(train_size, sizeof(double));
  double* Bixi_mean = new double[train_size];
  double* z = new double[train_size];
  double* b_tmp;

  if (solution_method == CG)
    b_tmp = (double*)calloc(subproblem_length, sizeof(double));

  //std::cout << "traning train size: " << train_size << std::endl;
  for (int i = 0; i < train_size; i++) {
    z[i] = classes.get(i) / processes;
    //if (rank == 0) std::cout << z[i] << std::endl;
  }

  double* u = new double[train_size];
  double* b = new double[train_size];
  double* error = new double[train_size];

  for (int i = 0; i < train_size; i++) {
    Bixi[i] = Bixi_mean[i] = u[i] = 0.0;
  }

  double* alpha = (double*)calloc(subproblem_length, sizeof(double));
  double* alpha_shared = (double*)calloc(subproblem_length, sizeof(double));
  double* alpha_total;

  if (rank == 0) {
    alpha_total = new double[grid_size];
  }

  double prediction_error = convergence_threshold * 100; // MSE
  double prev_prediction_error = -1.0;
  int info;
  sg::base::SGppStopwatch* iterStopwatch = new sg::base::SGppStopwatch();
  iterStopwatch->start();

  while (true) {
    admm_iteration += 1;

    //update alpha
    // calculate right handside
    for (int i = 0; i < train_size; i++) {
      b[i] = Bixi[i] - Bixi_mean[i] + z[i] - u[i];
    }

    // calculate new alphas
    if (solution_method == Cholesky) {
      cblas_dgemv(CblasRowMajor, CblasTrans, train_size, subproblem_length, rho, B_, subproblem_length, b, 1, 0.0, alpha, 1);
      this->dpotrs("U", subproblem_length, 1, system_matrix, subproblem_length, alpha, subproblem_length, info);
    } else {
      //TODO: check if b_tmp is really needed
      cblas_dgemv(CblasRowMajor, CblasTrans, train_size, subproblem_length, rho, B_, subproblem_length, b, 1, 0.0, b_tmp, 1);
      this->solve_cg(b_tmp, alpha, info);
    }

    //if (info != 0) throw new sg::base::solver_exception("solver failed");
    // calculate Bixi
    // correct
    cblas_dgemv(CblasRowMajor, CblasNoTrans, train_size, subproblem_length, 1.0, B_, subproblem_length, alpha, 1, 0.0, Bixi, 1);

    // calculate Bixi_mean
    MPI_Allreduce(Bixi, Bixi_mean, train_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (num_shared_functions > 0) {
      memcpy(alpha_shared, alpha, num_shared_functions * sizeof(double));
      //calculate shared part of the result
      cblas_dgemv(CblasRowMajor, CblasNoTrans, train_size, subproblem_length, 1.0, B_, subproblem_length, alpha_shared, 1, 0.0, Bixi_shared, 1);

      MPI_Allreduce(Bixi_shared, Bixi_shared_sum, train_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      cblas_dscal(train_size, 1.0 / processes, Bixi_shared_sum, 1 );

      // Dirk's formula: Bixi_mean = 1/num_tasks (\sum_{i=1}{num_tasks} Bixi - (num_tasks-1)/num_tasks \sum_{i=1}{num_tasks} Bixi_shared)
      cblas_daxpy(train_size, -((double)(processes - 1))/*/processes*/, Bixi_shared_sum, 1, Bixi_mean, 1);
    }

    //estimate new error
    cblas_dcopy(train_size, Bixi_mean, 1, error, 1);
    cblas_daxpy(train_size, -1.0, classes.getPointer(), 1, error, 1);
    prediction_error = cblas_dnrm2(train_size, error, 1);
    prediction_error = prediction_error * prediction_error / train_size;

    if (rank == 0 && (admm_iteration < 100 || admm_iteration % 100 == 0))
      std::cout << "iteration  " << admm_iteration << " new prediction error " << prediction_error << std::endl;

    if (prediction_error < convergence_threshold || prev_prediction_error == prediction_error ) //fabs(prev_prediction_error - prediction_error) < 0.00001 * convergence_threshold)
      break;

    // update z = 1/(N+rho)(y+rho*(Bixi_mean + u))
    // calculate the mean of Bixi
    cblas_dscal(train_size, 1.0 / processes, Bixi_mean, 1);
    // c = Bixi_mean + u
    cblas_daxpy(train_size, 1.0, Bixi_mean, 1, u, 1);
    cblas_dcopy(train_size, classes.getPointer(), 1, z, 1);
    cblas_daxpy(train_size, rho, u, 1, z, 1);
    cblas_dscal(train_size, 1.0 / (processes + rho), z, 1);

    // update u
    cblas_daxpy(train_size, -1.0, z, 1, u, 1);

    prev_prediction_error = prediction_error;
  }

  if (rank == 0)
    std::cout << "average iteration time " << (iterStopwatch->stop() / admm_iteration) << " seconds" << std::endl;

  int* rcounts, *displs;

  if (rank == 0) {
    rcounts = new int[processes];
    displs = new int[processes];
    int grid_points = this->grid_->getSize();
    int rank_subproblem_length;
    //int part_size = grid_points / num_tasks + ((grid_points%num_tasks) ? 1 : 0);
    int part_size = subproblem_length - num_shared_functions;

    for (int i = 0; i < processes; i++) {
      /*first_index = part_size * i;
      last_index = (num_shared_functions + part_size * (i+1) < grid_points) ? part_size * (i+1) - 1 : grid_points - 1;
      rank_subproblem_length = last_index - first_index + 1;*/


      int part_size = (grid_points - num_shared_functions) / num_tasks + (((grid_points - num_shared_functions) % num_tasks) ? 1 : 0);
      int fi = part_size * i + num_shared_functions;
      int li = (part_size * (i + 1) + num_shared_functions < grid_points) ? part_size * (i + 1) - 1 + num_shared_functions : grid_points - 1;
      rank_subproblem_length = li - fi + 1  + num_shared_functions;

      rcounts[i] = rank_subproblem_length;
      displs[i] = first_index;
    }
  }

  // collect surplusses
  MPI_Gatherv(alpha + num_shared_functions, subproblem_length - num_shared_functions, MPI_DOUBLE, alpha_total + num_shared_functions, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  memcpy(alpha_total, alpha, num_shared_functions);
  result.timeComplete_ =  myStopwatch->stop();

  if (rank == 0) {
    std::cout << "Final MSE: " << prediction_error << std::endl;
    std::cout << "Total number of ADMM operations: " << admm_iteration << std::endl;
    std::cout << std::endl;
  }

  delete [] Bixi, Bixi_mean, z, u, b, error;
  free(Bixi_shared);
  free(Bixi_shared_sum);
  delete myStopwatch;

  if (rank == 0) {
    alpha_ = new DataVector(alpha_total, grid_size);
    delete [] alpha_total, rcounts, displs;
  }

  if (solution_method == CG)
    free(b_tmp);

  free(alpha);
  free(alpha_shared);
  return result;

}

