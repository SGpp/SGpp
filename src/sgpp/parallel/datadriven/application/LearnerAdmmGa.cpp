/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "parallel/datadriven/application/LearnerBlockGaussSeidel.hpp"
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
#include <mpi.h>
//#include "omp.h"


#include "datadriven/algorithm/DMSystemMatrixExplicit.hpp"


using namespace sg::parallel;
using namespace sg::datadriven;
using namespace sg::base;


/**
 * More canonical ADMM implementation. Bixi_mean and z values used in the calculation
 * are "average" values.
 */
LearnerTiming LearnerBlockGaussSeidel::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
    const double lambda, double convergence_threshold, double scaling_factor) {
  LearnerTiming result;
  sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
  myStopwatch->start();

  // initializations
  this->createFactorizedSystem(trainDataset, lambda);
  int grid_size = this->grid_->getSize();

  if (scaling_factor == -1) scaling_factor = num_tasks;

  double* Bixi = new double[train_size]; //(double*)std::calloc(train_size, sizeof(double));
  double* Bixi_mean = new double[train_size];
  double* z = new double[train_size];
  double* b_tmp = (double*)calloc(subproblem_length, sizeof(double));

  for (int i = 0; i < train_size; i++) {
    // *** This part is different prom ADMM***
    z[i] = classes.get(i) / scaling_factor;
  }

  double* u = new double[train_size];
  double* b = new double[train_size];
  double* error = new double[train_size];

  for (int i = 0; i < train_size; i++) {
    Bixi[i] = Bixi_mean[i] = u[i] = 0.0;
  }

  double* alpha = (double*)calloc(subproblem_length, sizeof(double));
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
    // ***This part is different from ADMM***
    // calculate right handside
    if (rank > 0 || admm_iteration > 1) {
      MPI_Recv(Bixi_mean, train_size, MPI_DOUBLE, ((rank > 0) ? rank - 1 : num_tasks - 1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //cblas_daxpy(train_size, -1.0, Bixi, 1, Bixi_mean, 1);
    }

    for (int i = 0; i < train_size; i++) {
      b[i] = Bixi[i] - 1.0 / scaling_factor * Bixi_mean[i] + z[i] - u[i];
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

    if (admm_iteration > 1) cblas_daxpy(train_size, -1.0, Bixi, 1, Bixi_mean, 1);

    //if (info != 0) throw new sg::base::solver_exception("solver failed");
    // calculate Bixi
    cblas_dgemv(CblasRowMajor, CblasNoTrans, train_size, subproblem_length, 1.0, B_, subproblem_length, alpha, 1, 0.0, Bixi, 1);
    // calculate Bixi_mean
    //MPI_Allreduce(Bixi, Bixi_mean, train_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    cblas_daxpy(train_size, 1.0, Bixi, 1, Bixi_mean, 1);
    MPI::COMM_WORLD.Isend(Bixi_mean, train_size, MPI_DOUBLE, (rank + 1) % num_tasks, 0);

    //MPI_Allreduce(Bixi, Bixi_mean, train_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    //estimate new error
    cblas_dcopy(train_size, Bixi_mean, 1, error, 1);
    cblas_daxpy(train_size, -1.0, classes.getPointer(), 1, error, 1);
    prediction_error = cblas_dnrm2(train_size, error, 1);
    prediction_error = prediction_error * prediction_error / train_size;

    if (rank == 0 && (admm_iteration < 100 || admm_iteration % 100 == 0))
      std::cout << "iteration  " << admm_iteration << " new prediction error " << prediction_error << std::endl;

    if (prediction_error < convergence_threshold || prev_prediction_error == prediction_error)
      break;

    // update z
    //cblas_dscal(train_size, 1.0/processes, Bixi_mean, 1);
    // c = Bixi_mean + u
    cblas_daxpy(train_size, 1.0 / scaling_factor, Bixi_mean, 1, u, 1);
    cblas_dcopy(train_size, classes.getPointer(), 1, z, 1);
    cblas_daxpy(train_size, rho, u, 1, z, 1);
    cblas_dscal(train_size, 1.0 / (scaling_factor + rho), z, 1);

    // update u
    cblas_daxpy(train_size, -1.0, z, 1, u, 1);

    prev_prediction_error = prediction_error;
  }

  if (rank == 0)
    std::cout << "average iteration time " << (iterStopwatch->stop() / admm_iteration) << " seconds" << std::endl;

  int* rcounts, *displs;

  if (rank == 0) {
    rcounts = new int[num_tasks];
    displs = new int[num_tasks];
    int grid_points = this->grid_->getSize();
    int rank_subproblem_length;
    int part_size = grid_points / num_tasks + ((grid_points % num_tasks) ? 1 : 0);

    for (int i = 0; i < num_tasks; i++) {
      first_index = part_size * i;
      last_index = (part_size * (i + 1) < grid_points) ? part_size * (i + 1) - 1 : grid_points - 1;
      rank_subproblem_length = last_index - first_index + 1;
      rcounts[i] = rank_subproblem_length;
      displs[i] = first_index;
    }
  }

  // collect surplusses
  MPI_Gatherv(alpha, subproblem_length, MPI_DOUBLE, alpha_total, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  result.timeComplete_ =  myStopwatch->stop();

  if (rank == 0) {
    std::cout << "Final MSE: " << prediction_error << std::endl;
    std::cout << "Total number of ADMM operations: " << admm_iteration << std::endl;
    std::cout << std::endl;
  }

  delete [] Bixi, Bixi_mean, z, u, b, error;
  delete myStopwatch;

  if (rank == 0) {
    delete [] alpha_total, rcounts, displs;
  }

  free(b_tmp);
  free(alpha);
  return result;

}

