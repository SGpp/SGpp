/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "parallel/datadriven/application/LearnerAdmmCt.hpp"
#include "datadriven/algorithm/DMSystemMatrixExplicitFullgridBoundary.hpp"
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
#include "combigrid/combischeme/CombiS_CT.hpp"
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

extern "C"
{
  int dpotrf_(const char* uplo, const int* n, double* a, const int *
              lda, int* info);
  int dpotrs_(const char* uplo, const int* n, const int* nrhs, double* a,
              const int* lda, double* b, const int* ldb, int* info);
}

int LearnerAdmmCt::dpotrf(char* uplo, int n, double* a, int lda, int& info) {
  return dpotrf_(uplo, &n, a, &lda, &info);
}

int LearnerAdmmCt::dpotrs(const char* uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info) {
  return dpotrs_(uplo, &n, &nrhs, a, &lda, b, &ldb,  &info);
}

void LearnerAdmmCt::createFactorizedSystem(DataMatrix& trainDataset, double lambda) {
  double rho = 1.0;
  int dim = trainDataset.getNcols();

  train_size = trainDataset.getNrows();

  DMSystemMatrixExplicitFullgridBoundary explicit_matrix(trainDataset, lambda, level_vector, dim);
  grid_size = explicit_matrix.get_size();
  explicit_matrix.compileMatrix(&system_matrix, &B_);

  //perform Cholesky factorization
  int info;
  dpotrf(const_cast<char*>("U"), grid_size, system_matrix, grid_size, info);

  if (info != 0) {
    std::cerr << "Cholesky factorization failed" << std::endl;
    throw new std::exception();
  }
}


LearnerAdmmCt::LearnerAdmmCt(const sg::base::RegularGridConfiguration& GridConfig, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose )
  : LearnerAdmm(regularization, isRegression, isVerbose), rho(1.0), admm_iteration(0) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  B_ = NULL;
  system_matrix = NULL;

  InitializeGrid(GridConfig);

}

void LearnerAdmmCt::InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig) {
  combigrid::S_CT ct(GridConfig.dim_, GridConfig.level_);
  level_vector = ct.getLevel(rank);
  comb_coeff = ct.getCoef(rank);
}

/**
 * Destructor
 */
LearnerAdmmCt::~LearnerAdmmCt() {
  if (system_matrix != NULL)
    delete[] system_matrix;

  if (B_ != NULL)
    delete[] B_;
}

LearnerTiming LearnerAdmmCt::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
                                   const double lambda, double convergence_threshold) {

  LearnerTiming result;
  sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
  myStopwatch->start();

  // initializations
  this->createFactorizedSystem(trainDataset, lambda);
  int processes = num_tasks;
  double* Bixi = new double[train_size]; //(double*)std::calloc(train_size, sizeof(double));
  double* Bixi_mean = new double[train_size];
  double* z = new double[train_size];

  for (int i = 0; i < train_size; i++) {
    z[i] = classes.get(i) / processes;
  }

  double* u = new double[train_size];
  double* b = new double[train_size];
  double* error = new double[train_size];

  for (int i = 0; i < train_size; i++) {
    Bixi[i] = Bixi_mean[i] = u[i] = 0.0;
  }

  double* alpha = new double[grid_size];
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
    //TODO: it should be implemented more efficiently
    for (int i = 0; i < train_size; i++) {
      b[i] = comb_coeff * Bixi[i] - Bixi_mean[i] + z[i] - u[i];
    }


    cblas_dgemv(CblasRowMajor, CblasTrans, train_size, grid_size, rho, B_, grid_size, b, 1, 0.0, alpha, 1);

    // calculate new alphas
    this->dpotrs("U", grid_size, 1, system_matrix, grid_size, alpha, grid_size, info);

    if (info != 0)
      throw new sg::base::solver_exception("Cholesky solver failed");

    // calculate Bixi
    cblas_dgemv(CblasRowMajor, CblasNoTrans, train_size, grid_size, rho, B_, grid_size, alpha, 1, 0.0, Bixi, 1);
    //multiply Bixi by combination coefficient (to produce combigrid solution as suggested by Matthias)
    cblas_dscal(train_size, comb_coeff, Bixi, 1);
    // calculate Bixi_mean
    MPI_Allreduce(Bixi, Bixi_mean, train_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //estimate new error
    cblas_dcopy(train_size, Bixi_mean, 1, error, 1);
    cblas_daxpy(train_size, -1.0, classes.getPointer(), 1, error, 1);
    prediction_error = cblas_dnrm2(train_size, error, 1);
    prediction_error = prediction_error * prediction_error / train_size;

    if (rank == 0 && (admm_iteration < 100 || admm_iteration % 100 == 0))
      std::cout << "iteration  " << admm_iteration << " new prediction error " << prediction_error << std::endl;

    if (prediction_error < convergence_threshold || fabs(prev_prediction_error - prediction_error) < 0.001 * convergence_threshold)
      break;

    // update z
    cblas_dscal(train_size, 1.0 / processes, Bixi_mean, 1);
    //cblas_dscal(train_size, 1.0, Bixi_mean, 1);
    // c = Bixi_mean + u
    cblas_daxpy(train_size, 1.0, Bixi_mean, 1, u, 1);
    //cblas_dcopy(train_size, classes.getPointer(), 1, z, 1);
    cblas_daxpy(train_size, rho, u, 1, z, 1);
    cblas_dscal(train_size, 1.0 / (processes + rho), z, 1);

    // update u
    cblas_daxpy(train_size, -1.0, z, 1, u, 1);

    prev_prediction_error = prediction_error;
  }

  if (rank == 0)
    std::cout << "average iteration time " << (iterStopwatch->stop() / admm_iteration) << " seconds" << std::endl;

  /*
   //collect results
  int *rcounts, *displs;
  if (rank == 0)
  {
      rcounts = new int[processes];
      displs = new int[processes];
      int grid_size = this->grid_->getSize();
      int rank_subproblem_length;
      int part_size = grid_size / num_tasks + ((grid_size%num_tasks) ? 1 : 0);
      for (int i = 0; i < processes; i++)
      {
          start_index = part_size * i;
          end_index = (part_size * (i+1) < grid_size) ? part_size * (i+1) - 1 : grid_size - 1;
          rank_subproblem_length = end_index - start_index + 1;
          rcounts[i] = rank_subproblem_length;
          displs[i] = start_index;
      }
  }

  // collect surplusses
  MPI_Gatherv(alpha, grid_size, MPI_DOUBLE, alpha_total, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
  result.timeComplete_ =  myStopwatch->stop();

  if (rank == 0) {
    std::cout << "Final MSE: " << prediction_error << std::endl;
    std::cout << "Total number of ADMM operations: " << admm_iteration << std::endl;
    std::cout << std::endl;
  }

  delete [] Bixi, Bixi_mean, z, u, b, error, alpha;
  delete myStopwatch;

  if (rank == 0) {
    delete [] alpha_total/*, rcounts, displs*/;
  }

  return result;

}

