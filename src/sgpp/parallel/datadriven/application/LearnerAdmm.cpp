/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "parallel/datadriven/application/LearnerAdmm.hpp"
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

extern "C"
{
  int dpotrf_(const char* uplo, const int* n, double* a, const int *
              lda, int* info);
  int dpotrs_(const char* uplo, const int* n, const int* nrhs, double* a,
              const int* lda, double* b, const int* ldb, int* info);
}

int LearnerAdmm::dpotrf(char* uplo, int n, double* a, int lda, int& info) {
  return dpotrf_(uplo, &n, a, &lda, &info);
}

int LearnerAdmm::dpotrs(const char* uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info) {
  return dpotrs_(uplo, &n, &nrhs, a, &lda, b, &ldb,  &info);
}


void LearnerAdmm::InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig) {
  if (GridConfig.type_ == sg::base::LinearTrapezoidBoundary) {
    grid_ = new sg::base::LinearTrapezoidBoundaryGrid(GridConfig.dim_);
  } else if (GridConfig.type_ == sg::base::ModLinear) {
    grid_ = new sg::base::ModLinearGrid(GridConfig.dim_);
  } else if (GridConfig.type_ == sg::base::Linear) {
    grid_ = new sg::base::LinearGrid(GridConfig.dim_);
  } else {
    grid_ = NULL;
    throw base::application_exception("LearnerADMM::InitializeGrid: An unsupported grid type was chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  sg::base::GridGenerator* myGenerator = grid_->createGridGenerator();
  myGenerator->regular(GridConfig.level_);

  delete myGenerator;
}


void LearnerAdmm::createFactorizedSystem(DataMatrix& trainDataset, double lambda) {
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

  if (isVerbose_)
    std::cout << "rank " << rank << " subproblem size " << subproblem_length << std::endl;

  if (isVerbose_ && rank == 0)
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
  if (isVerbose_ && rank == 0)
    std::cout << "lambda value: " << lambda << std::endl;

  cblas_dsyrk(CblasRowMajor, CblasLower,
              CblasTrans, cols, rows, rho, B_,
              cols, lambda, system_matrix, cols);

  if (isVerbose_)
    std::cout << "rank " << rank << " building matrix directly in " << myStopwatch->stop() << " seconds" << std::endl;

  //Cholesky factorization
  myStopwatch->start();

  int info;

  if (solution_method == Cholesky) {
    this->dpotrf("U", subproblem_length, system_matrix, subproblem_length, info);

    if (isVerbose_)
      std::cout << "rank " << rank << " performing cholesky factorization in " << myStopwatch->stop() << " seconds" << std::endl;
  }

  delete myStopwatch;
}

void LearnerAdmm::set_rho(double rho) {
  this->rho = rho;
}
double LearnerAdmm::get_rho() {
  return this->rho;
}

LearnerAdmm::LearnerAdmm(const sg::base::RegularGridConfiguration& GridConfig, sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose, SolutionType method )
  : Learner(regularization, isRegression, isVerbose), rho(1.0), admm_iteration(0), solution_method(method) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  this->InitializeGrid(GridConfig);
  init();
}

LearnerAdmm::LearnerAdmm(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose, SolutionType method )
  : Learner(regularization, isRegression, isVerbose), rho(1.0), admm_iteration(0), solution_method(method) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  init();
}

LearnerAdmm::LearnerAdmm(std::string tGridFilename, std::string tAlphaFilename, sg::datadriven::LearnerRegularizationType& regularization,
                         const bool isRegression, const bool isVerbose, SolutionType method)
  : Learner(tGridFilename, tAlphaFilename, regularization, isRegression, isVerbose), rho(1.0), admm_iteration(0), solution_method(method) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  init();
}

void LearnerAdmm::init() {
  B_ = NULL;
  C_ = NULL;
  system_matrix = NULL;
  alpha_ = NULL;
  cg_residual_norm_threshold = 0.00001;
}

/**
 * Destructor
 */
LearnerAdmm::~LearnerAdmm() {
  /*if (C_ != NULL)
      delete C_;*/

  if (system_matrix != NULL)
    delete[] system_matrix;

  if (B_ != NULL)
    delete[] B_;

  /*if (grid_ != NULL)
      delete grid_;

  if (alpha_ != NULL)
      delete alpha_;*/
}

LearnerTiming LearnerAdmm::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
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
    //if (rank == 0) std::cout << z[i] << std::endl;
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

    //estimate new error
    cblas_dcopy(train_size, Bixi_mean, 1, error, 1);
    cblas_daxpy(train_size, -1.0, classes.getPointer(), 1, error, 1);
    prediction_error = cblas_dnrm2(train_size, error, 1);
    prediction_error = prediction_error * prediction_error / train_size;

    if (isVerbose_ && rank == 0 && (admm_iteration < 100 || admm_iteration % 100 == 0))
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

  if (isVerbose_ && rank == 0)
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

// do a several iterations of CG on Ax = b
void LearnerAdmm::solve_cg(double* b, double* x, int& info) {
  int imax = subproblem_length;

  //double* t = (double *)calloc(subproblem_length, sizeof(double));
  // compute r = b - Ax
  double* r = (double*)malloc(subproblem_length * sizeof(double));
  cblas_dcopy(subproblem_length, b, 1, r, 1);
  cblas_dsymv(CblasRowMajor, CblasLower, subproblem_length,  -1.0, system_matrix, subproblem_length, x, 1, 1.0, r, 1);

  bool converged = false;
  double rho_new, rho_old, residual_norm, beta, alpha;
  double* p = (double*)malloc(subproblem_length * sizeof(double));
  double* q = (double*)malloc(subproblem_length * sizeof(double));

  // for i = 1, 2, ...
  for (int i = 1; !converged && i < imax; i++) {
    // \rho\{i-1} = r^{(i-1)}^Tr^{(i-1)}
    rho_new = cblas_dnrm2(subproblem_length, r, 1);
    rho_new *= rho_new;

    // if i = 1
    if (i == 1) {
      // p^{(1)} = r^{(0)}
      cblas_dcopy(subproblem_length, r, 1, p, 1);
    }
    // else
    else {
      //\beta_{i-1} = \rho_{i-1}/\rho_{i-2}
      beta = rho_new / rho_old;
      //p^{(i)} = r^{(i-1)} + \beta_{i-1} p^{(i-1)}
      cblas_dscal(subproblem_length, beta, p, 1);
      cblas_daxpy(subproblem_length, 1.0, r, 1, p, 1);
    }

    // end if
    // q^{(i)} = A p^{(i)}
    cblas_dsymv(CblasRowMajor, CblasLower,  subproblem_length, 1.0, system_matrix, subproblem_length, p, 1, 0.0, q, 1);

    // alpha_i = \rho_{i-1}/(p^{(i)}^T q^{(i)}
    alpha = cblas_ddot(subproblem_length, p, 1, q, 1);
    alpha = rho_new / alpha;

    // x^{(i)} = x^{(i-1)} + \alpha_i p^{(i)}
    cblas_daxpy(subproblem_length, alpha, p, 1, x, 1);
    // r^{(i)} = r^{(i-1)} - \alpha_i q^{(i)}
    cblas_daxpy(subproblem_length, -alpha, q, 1, r, 1);

    //check convergence; continue if necessary
    residual_norm = cblas_dnrm2(subproblem_length, r, 1);

    //std::cout << "residual " << residual_norm << std::endl;
    if (isnan(residual_norm) ) {
      std::cerr << "residual_norm is not a number" << std::endl;
      break;
    }

    if (residual_norm < cg_residual_norm_threshold) {
      converged = true;
    }

    rho_old = rho_new;
  }

  info = converged ? 0 : -1;
  //end for
  free(r);
  free(p);
  free(q);
}
