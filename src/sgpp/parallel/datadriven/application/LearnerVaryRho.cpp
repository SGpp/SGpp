/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "parallel/datadriven/application/LearnerVaryRho.hpp"
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
#include <math.h>
//#include "omp.h"

#include "datadriven/algorithm/DMSystemMatrixExplicit.hpp"


using namespace sg::parallel;
using namespace sg::datadriven;
using namespace sg::base;

void LearnerVaryRho::createFactorizedSystem(DataMatrix& trainDataset, double lambda, bool verbose) {
  //double rho = 1.0;
  int grid_points = this->grid_->getSize();

  train_size = trainDataset.getNrows();
  int rows = train_size;

  // ceiling
  int part_size = grid_points / num_tasks + ((grid_points % num_tasks) ? 1 : 0);
  first_index = part_size * rank;
  last_index = (part_size * (rank + 1) < grid_points) ? part_size * (rank + 1) - 1 : grid_points - 1;
  subproblem_length = last_index - first_index + 1;
  int cols = subproblem_length;

  if (verbose) std::cout << "rank " << rank << " subproblem size " << subproblem_length << std::endl;

  if (rank == 0 && verbose)
    std::cout << "building matrix" << std::endl;

  //Assemble matrix
  sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
  myStopwatch->start();

  DMSystemMatrixExplicit dmsystem(trainDataset, lambda);

  if (B_ != NULL)
    delete B_;

  // correct
  B_ = dmsystem.compileMatrix(this->grid_->getStorage(), first_index, last_index);

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
  if (rank == 0 && verbose)
    std::cout << "lambda value: " << lambda << std::endl;

  cblas_dsyrk(CblasRowMajor, CblasLower,
              CblasTrans, cols, rows, rho, B_,
              cols, lambda, system_matrix, cols);

  if (verbose) std::cout << "rank " << rank << " building matrix directly in " << myStopwatch->stop() << " seconds" << std::endl;

  //Cholesky factorization
  myStopwatch->start();

  int info;

  if (solution_method == Cholesky) {
    this->dpotrf("U", subproblem_length, system_matrix, subproblem_length, info);

    if (verbose) std::cout << "rank " << rank << " performing cholesky factorization in " << myStopwatch->stop() << " seconds" << std::endl;
  }

  delete myStopwatch;
}

LearnerTiming LearnerVaryRho::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
                                    const double lambda, double convergence_threshold) {
  LearnerTiming result;
  sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
  myStopwatch->start();

  // initializations
  this->createFactorizedSystem(trainDataset, lambda);
  int grid_size = this->grid_->getSize();
  int processes = num_tasks;
  double* Bixi = new double[train_size]; //(double*)std::calloc(train_size, sizeof(double));
  double* Bixi_mean = new double[train_size];
  double* z = new double[train_size];
  double* b_tmp;

  if (solution_method == CG)
    b_tmp = (double*)calloc(subproblem_length, sizeof(double));

  //std::cout << "traning train size: " << train_size << std::endl;
  for (int i = 0; i < train_size; i++) {
    z[i] = classes.get(i) / processes;
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

  double tau_incr = 2;
  double tau_decr = 2;
  double mu = 20;
  double* residual_primal, *residual_dual, *z_old;
  double rho_old, s_norm, r_norm, s_norm_total, r_norm_total;
  s_norm_total = r_norm_total = 0.0;
  residual_primal = (double*)malloc(train_size * sizeof(double));
  residual_dual = (double*)calloc(subproblem_length, sizeof(double));
  z_old = (double*)malloc(train_size * sizeof(double));

  sg::base::SGppStopwatch* iterStopwatch = new sg::base::SGppStopwatch();
  iterStopwatch->start();

  while (true) {
    rho_old = rho;

    if (r_norm > mu * s_norm )
      rho = tau_incr * rho;
    else if (s_norm > mu * r_norm || (admm_iteration > 5 && prev_prediction_error - prev_prediction_error < 0.5))
      rho = rho / tau_decr;

    /*else
        rho = rho;*/

    if (rank == 0) std::cout << "new rho " << rho << std::endl;

    if (rho_old != rho) createFactorizedSystem(trainDataset, lambda);


    admm_iteration += 1;

    //update alpha
    // calculate right handside
    for (int i = 0; i < train_size; i++) {
      b[i] = Bixi[i] - Bixi_mean[i] + z[i] - u[i];
    }

    // calculate new alphas
    if (solution_method == Cholesky) {
      cblas_dgemv(CblasRowMajor, CblasTrans, train_size, subproblem_length, rho, B_, subproblem_length, b, 1, 0.0, alpha, 1);
      //if(rank == 0) std::cout << "b: " << cblas_dnrm2(subproblem_length, alpha, 1) << std::endl;

      this->dpotrs("U", subproblem_length, 1, system_matrix, subproblem_length, alpha, subproblem_length, info);
      //if (rank == 0) std::cout << "alpha norm: " << cblas_dnrm2(subproblem_length, alpha, 1) << std::endl;

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
    cblas_dcopy(train_size, z, 1, z_old, 1);
    cblas_dscal(train_size, 1.0 / processes, Bixi_mean, 1);
    // c = Bixi_mean + u
    cblas_daxpy(train_size, 1.0, Bixi_mean, 1, u, 1);
    cblas_dcopy(train_size, classes.getPointer(), 1, z, 1);
    cblas_daxpy(train_size, rho, u, 1, z, 1);
    cblas_dscal(train_size, 1.0 / (processes + rho), z, 1);

    // update u
    cblas_daxpy(train_size, -1.0, z, 1, u, 1);

    // udpate primal residual r = (Bixi_mean - z)
    cblas_dcopy(train_size, Bixi_mean, 1, residual_primal, 1);
    cblas_daxpy(train_size, -1, z  , 1, residual_primal, 1);
    r_norm = cblas_dnrm2(train_size, residual_primal, 1);

    /*r_norm = r_norm * r_norm;
    MPI_Allreduce(&r_norm, &r_norm_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r_norm_total = sqrt(r_norm_total);*/
    if (rank == 0) std::cout << "new r_norm: " << r_norm << std::endl;


    //update dual residual s = -rho Bi^T I (z_old - z_new)
    cblas_daxpy(train_size, -1.0, z, 1, z_old, 1);
    cblas_dgemv(CblasRowMajor, CblasTrans, train_size, subproblem_length, -rho, B_, subproblem_length, z_old, 1, 0, residual_dual, 1);
    s_norm = cblas_dnrm2(subproblem_length, residual_dual, 1);

    /*s_norm = s_norm * s_norm;
    MPI_Allreduce(&s_norm, &s_norm_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    s_norm_total = sqrt(s_norm_total);*/
    if (rank == 0) std::cout << "new s_norm: " << s_norm << std::endl;

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
    int part_size = grid_points / num_tasks + ((grid_points % num_tasks) ? 1 : 0);

    for (int i = 0; i < processes; i++) {
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
    alpha_ = new DataVector(alpha_total, grid_size);
    delete [] alpha_total, rcounts, displs;
  }

  if (solution_method == CG)
    free(b_tmp);

  free(alpha);
  return result;

}

