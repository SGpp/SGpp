// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt() : DBMatOffline() {
  this->q_ortho_matrix_ = sgpp::base::DataMatrix(1, 1, 1.0);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(1, 1);
  // Deprecated
  // dim_a = 0, indirectly tells the online object if build() or decompose() were performed
}

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const std::string& fileName)
    : DBMatOffline(fileName) {
  // Read grid size from header (number of rows in lhsMatrix)
  std::ifstream filestream(fileName, std::istream::in);
  // Read configuration
  if (!filestream) {
    throw sgpp::base::algorithm_exception("Failed to open File");
  }
  std::string str;
  std::getline(filestream, str);
  filestream.close();

  std::vector<std::string> tokens;
  StringTokenizer::tokenize(str, ",", tokens);

  auto size = std::stoi(tokens[0]);
  std::cout << "Grid size " << size << std::endl;
  // grid already initialized in super constructor

  this->lhsMatrix = sgpp::base::DataMatrix(size, size);
  this->q_ortho_matrix_ = sgpp::base::DataMatrix(size, size);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(size, size);
#ifdef USE_GSL
  gsl_matrix_view lhs_view = gsl_matrix_view_array(this->lhsMatrix.getPointer(), size, size);
  gsl_matrix_view q_view = gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), size, size);
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(this->t_tridiag_inv_matrix_.getPointer(), size, size);

  FILE* file = fopen(fileName.c_str(), "rb");
  if (!file) {
    throw sgpp::base::algorithm_exception{"Failed to open File"};
  }

  // seek end of first line
  char c = 0;
  while (c != '\n') {
    c = static_cast<char>(fgetc(file));
  }

  std::cout << "Fread init" << std::endl;
  gsl_matrix_fread(file, &lhs_view.matrix);
  gsl_matrix_fread(file, &q_view.matrix);
  gsl_matrix_fread(file, &t_inv_view.matrix);
  std::cout << "Fread done" << std::endl;

  fclose(file);

  this->isConstructed = true;
  this->isDecomposed = true;
#else
  throw sgpp::base::algorithm_exception("USE_GSL has to be set");
#endif /* USE_GSL */
}

DBMatOffline* DBMatOfflineOrthoAdapt::clone() { return new DBMatOfflineOrthoAdapt{*this}; }

bool DBMatOfflineOrthoAdapt::isRefineable() { return true; }

void DBMatOfflineOrthoAdapt::buildMatrix(Grid* grid,
                                         RegularizationConfiguration& regularizationConfig) {
  DBMatOffline::buildMatrix(grid, regularizationConfig);
  size_t dim_a = grid->getStorage().getSize();

  this->q_ortho_matrix_.resizeQuadratic(dim_a);
  this->t_tridiag_inv_matrix_.resizeQuadratic(dim_a);
  // isConstructed = true, set by parent call
}

void DBMatOfflineOrthoAdapt::decomposeMatrix(
    RegularizationConfiguration& regularizationConfig,
    DensityEstimationConfiguration& densityEstimationConfig) {
#ifdef USE_GSL
  if (!isConstructed) {
    throw sgpp::base::algorithm_exception(
        "in DBMatOfflineOrthoAdapt::decomposeMatrix: \nmatrix not built yet.");
  }
  size_t dim_a = lhsMatrix.getNrows();
  if (dim_a <= 1) {
    this->q_ortho_matrix_.set(0, 0, 1.0);
    this->t_tridiag_inv_matrix_.set(0, 0, 1.0 / this->lhsMatrix.get(0, 0));
    isDecomposed = true;
    return;
  }
  // allocating subdiagonal and diagonal vectors of T
  sgpp::base::DataVector diag(dim_a);
  sgpp::base::DataVector subdiag(dim_a - 1);

  // decomposing: lhs = Q * T * Q^t
  this->hessenberg_decomposition(diag, subdiag);

  // adding configuration parameter lambda to diag before inverting T
  for (size_t i = 0; i < dim_a; i++) {
    diag.set(i, diag.get(i) + regularizationConfig.lambda_);
  }

  // inverting T+lambda*I, by solving L*R*x_i = e_i, for every i-th column x_i of T_inv
  this->invert_symmetric_tridiag(diag, subdiag);

  // decomposed matrix now consists of Q, T^-1, with: (lhs+lambda*I)^-1 = Q * T^-1 * Q^t
  this->isDecomposed = true;
#else
  throw base::not_implemented_exception("built without GSL");
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::decomposeMatrixParallel(
    RegularizationConfiguration& regularizationConfig,
    DensityEstimationConfiguration& densityEstimationConfig,
    std::shared_ptr<BlacsProcessGrid> processGrid, const ParallelConfiguration& parallelConfig) {
#ifdef USE_SCALAPACK
  if (!isConstructed) {
    throw sgpp::base::algorithm_exception(
        "in DBMatOfflineOrthoAdapt::decomposeMatrix: \nmatrix not built yet.");
  }

  // syncing lhs distributed matrix
  this->lhsDistributed = DataMatrixDistributed::fromSharedData(
      this->lhsMatrix.getPointer(), processGrid, lhsMatrix.getNrows(), lhsMatrix.getNcols(),
      parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  size_t dim_a = lhsMatrix.getNrows();
  if (dim_a <= 1) {
    this->q_ortho_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
        this->q_ortho_matrix_.getPointer(), processGrid, this->q_ortho_matrix_.getNrows(),
        this->q_ortho_matrix_.getNcols(), parallelConfig.rowBlockSize_,
        parallelConfig.columnBlockSize_);
    this->q_ortho_matrix_distributed_.set(0, 0, 1.0);
    this->t_tridiag_inv_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
        this->t_tridiag_inv_matrix_.getPointer(), processGrid,
        this->t_tridiag_inv_matrix_.getNrows(), this->t_tridiag_inv_matrix_.getNcols(),
        parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);
    this->t_tridiag_inv_matrix_distributed_.set(0, 0, 1.0 / this->lhsMatrix.get(0, 0));
    isDecomposed = true;
    return;
  }

  // allocating subdiagonal and diagonal vectors of T
  sgpp::base::DataVector d(dim_a);
  sgpp::base::DataVector sd(dim_a - 1);
  sgpp::base::DataVector tau(dim_a);

  for (size_t i = 0; i < tau.size(); i++) {
    tau.set(i, 1.1);
  }
  int lwork = -1;
  double* work = new double[1];
  int info;
  // ask pdsytrd_ how much workspace it needs
  pdsytrd_("L", dim_a, this->lhsDistributed.getLocalPointer(), 1, 1,
           this->lhsDistributed.getDescriptor(), d.getPointer(), sd.getPointer(), tau.getPointer(),
           work, lwork, info);
  lwork = work[0];
  work = new double[lwork];
  // pdsytrd_ amounts to hessenberg_decomposition of non-parallel version
  pdsytrd_("L", dim_a, this->lhsDistributed.getLocalPointer(), 1, 1,
           this->lhsDistributed.getDescriptor(), d.getPointer(), sd.getPointer(), tau.getPointer(),
           work, lwork, info);
  free(work);

  // inverting middle matrix with added lambda: T <~ (T + lambda*I)^-1
  for (size_t i = 0; i < dim_a; i++) {
    d.set(i, d.get(i) + regularizationConfig.lambda_);
  }

  this->invert_symmetric_tridiag(d, sd);
  this->t_tridiag_inv_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
      this->t_tridiag_inv_matrix_.getPointer(), processGrid, this->t_tridiag_inv_matrix_.getNrows(),
      this->t_tridiag_inv_matrix_.getNcols(), parallelConfig.rowBlockSize_,
      parallelConfig.columnBlockSize_);

  // obtaining Q
  this->q_ortho_matrix_.setAll(0.0);
  for (size_t i = 0; i < dim_a; i++) {
    this->q_ortho_matrix_.set(i, i, 1.0);
  }
  this->q_ortho_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
      this->q_ortho_matrix_.getPointer(), processGrid, this->q_ortho_matrix_.getNrows(),
      this->q_ortho_matrix_.getNcols(), parallelConfig.rowBlockSize_,
      parallelConfig.columnBlockSize_);
  // ask pdormtr_ how much workspace it needs
  double* work2 = new double[1];
  lwork = -1;
  pdormtr_("L", "L", "N", dim_a, dim_a, this->lhsDistributed.getLocalPointer(), 1, 1,
           this->lhsDistributed.getDescriptor(), tau.getPointer(),
           this->q_ortho_matrix_distributed_.getLocalPointer(), 1, 1,
           this->q_ortho_matrix_distributed_.getDescriptor(), work2, lwork, info);
  lwork = work2[0];
  work2 = new double[lwork];
  // pdormtr_ is used on Q (=identity) to obtain Q by overwriting Id*Q into it
  pdormtr_("L", "L", "N", dim_a, dim_a, this->lhsDistributed.getLocalPointer(), 1, 1,
           this->lhsDistributed.getDescriptor(), tau.getPointer(),
           this->q_ortho_matrix_distributed_.getLocalPointer(), 1, 1,
           this->q_ortho_matrix_distributed_.getDescriptor(), work2, lwork, info);
  free(work2);
  // sync non-distributed matrices
  q_ortho_matrix_distributed_.toLocalDataMatrix(q_ortho_matrix_);
  t_tridiag_inv_matrix_distributed_.toLocalDataMatrix(t_tridiag_inv_matrix_);

  this->isDecomposed = true;
#endif /* USE_SCALAPACK */
}

void DBMatOfflineOrthoAdapt::hessenberg_decomposition(sgpp::base::DataVector& diag,
                                                      sgpp::base::DataVector& subdiag) {
#ifdef USE_GSL
  size_t dim_a = lhsMatrix.getNrows();

  gsl_vector* tau = gsl_vector_alloc(dim_a - 1);
  gsl_matrix_view gsl_lhs = gsl_matrix_view_array(lhsMatrix.getPointer(), dim_a, dim_a);
  gsl_matrix_view gsl_q = gsl_matrix_view_array(q_ortho_matrix_.getPointer(), dim_a, dim_a);
  gsl_vector_view gsl_diag = gsl_vector_view_array(diag.getPointer(), dim_a);
  gsl_vector_view gsl_subdiag = gsl_vector_view_array(subdiag.getPointer(), dim_a - 1);

  // does the decomposition
  gsl_linalg_symmtd_decomp(&gsl_lhs.matrix, tau);

  // unpacks information out of matrix to explicitly create Q, and T
  gsl_linalg_symmtd_unpack(&gsl_lhs.matrix, tau, &gsl_q.matrix, &gsl_diag.vector,
                           &gsl_subdiag.vector);

  gsl_vector_free(tau);
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::invert_symmetric_tridiag(sgpp::base::DataVector& diag,
                                                      sgpp::base::DataVector& subdiag) {
#ifdef USE_GSL
  size_t dim_a = lhsMatrix.getNrows();
  if (dim_a <= 1) {
    return;
  }
  gsl_vector* e = gsl_vector_calloc(diag.getSize());  // calloc sets all values to zero
  gsl_vector* x = gsl_vector_alloc(diag.getSize());   // target of solving

  gsl_vector_view gsl_diag = gsl_vector_view_array(diag.getPointer(), dim_a);
  gsl_vector_view gsl_subdiag = gsl_vector_view_array(subdiag.getPointer(), dim_a - 1);

  // loops columns of T_inv
  for (size_t k = 0; k < t_tridiag_inv_matrix_.getNcols(); k++) {
    e->data[k] = 1;
    gsl_linalg_solve_symm_tridiag(&gsl_diag.vector, &gsl_subdiag.vector, e, x);
    for (size_t i = 0; i < t_tridiag_inv_matrix_.getNrows(); i++) {
      this->t_tridiag_inv_matrix_.set(k, i, x->data[i]);
    }
    e->data[k] = 0;
  }

  gsl_vector_free(e);
  gsl_vector_free(x);
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::store(const std::string& fileName) {
#ifdef USE_GSL
  DBMatOffline::store(fileName);

  FILE* outCFile = fopen(fileName.c_str(), "ab");
  if (!outCFile) {
    throw sgpp::base::algorithm_exception{"cannot open file for writing"};
  }

  auto dim_a = getGridSize();
  // store q_ortho_matrix_
  gsl_matrix_view q_view = gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), dim_a, dim_a);
  gsl_matrix_fwrite(outCFile, &q_view.matrix);

  // store t_inv_tridiag_
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(this->t_tridiag_inv_matrix_.getPointer(), dim_a, dim_a);
  gsl_matrix_fwrite(outCFile, &t_inv_view.matrix);

  fclose(outCFile);
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::syncDistributedDecomposition(
    std::shared_ptr<BlacsProcessGrid> processGrid, const ParallelConfiguration& parallelConfig) {
#ifdef USE_SCALAPACK
  if (!isDecomposed) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOfflineOrthoAdapt::syncDistributedDecomposition\nCan't sync, because lhsMatrix "
        "was not decomposed yet");
  }
  q_ortho_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
      q_ortho_matrix_.data(), processGrid, q_ortho_matrix_.getNrows(), q_ortho_matrix_.getNcols(),
      parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  t_tridiag_inv_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
      t_tridiag_inv_matrix_.data(), processGrid, t_tridiag_inv_matrix_.getNrows(),
      t_tridiag_inv_matrix_.getNcols(), parallelConfig.rowBlockSize_,
      parallelConfig.columnBlockSize_);
#endif
  // no action needed without scalapack
}

void DBMatOfflineOrthoAdapt::compute_inverse() {
#ifdef USE_GSL
  if (!isDecomposed) {
    throw sgpp::base::algorithm_exception(
        "in DBMatOfflineOrthoAdapt::compute_inverse:\noffline matrix not decomposed yet.\n");
  }

  // initialize lhsInverse
  this->lhsInverse = DataMatrix(this->lhsMatrix.getNrows(), this->lhsMatrix.getNcols());

  gsl_matrix_view inv_view = gsl_matrix_view_array(this->lhsInverse.getPointer(),
                                                   lhsInverse.getNrows(), lhsInverse.getNcols());

  gsl_matrix_view t_inv_view = gsl_matrix_view_array(
      this->getTinv().getPointer(), this->getTinv().getNrows(), this->getTinv().getNcols());

  gsl_matrix_view q_view = gsl_matrix_view_array(this->getQ().getPointer(), this->getQ().getNrows(),
                                                 this->getQ().getNcols());

  gsl_matrix* QT = gsl_matrix_alloc(lhsInverse.getNrows(), lhsInverse.getNcols());

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &q_view.matrix, &t_inv_view.matrix, 0.0, QT);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, QT, &q_view.matrix, 0.0, &inv_view.matrix);

  gsl_matrix_free(QT);
#else
  throw sgpp::base::algorithm_exception("build without GSL");
#endif /*USE_GSL*/
}

void DBMatOfflineOrthoAdapt::compute_inverse_parallel(std::shared_ptr<BlacsProcessGrid> processGrid,
                                                      const ParallelConfiguration& parallelConfig) {
#ifdef USE_SCALAPACK
  if (!isDecomposed) {
    throw sgpp::base::algorithm_exception(
        "in DBMatOfflineOrthoAdapt::compute_inverse_parallel:\noffline matrix not decomposed "
        "yet.\n");
  }
  size_t dim_a = this->lhsDistributed.getGlobalRows();

  if (dim_a <= 1) {
  }

  DataMatrixDistributed* QT = new DataMatrixDistributed(
      processGrid, dim_a, dim_a, parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);
  DataMatrixDistributed* INV = new DataMatrixDistributed(
      processGrid, dim_a, dim_a, parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  // Note: Because "L", a.k.a. lower triangular storage mode, was chosen on decomposition with
  // pdsytrd_, the Q is altered to Q^t, i.e. A = Q^t * T * Q, and also A^-1 = Q^t * T^-1 * Q
  DataMatrixDistributed::mult(this->q_ortho_matrix_distributed_,
                              this->t_tridiag_inv_matrix_distributed_, *QT, true, false);
  DataMatrixDistributed::mult(*QT, this->q_ortho_matrix_distributed_, *INV, false, false);

  // writing computed inverse into member and syncing both distri and non-distri matrices
  this->lhsInverse = DataMatrix(this->lhsMatrix.getNrows(), this->lhsMatrix.getNcols());
  INV->toLocalDataMatrix(this->lhsInverse);
  this->lhsDistributedInverse = DataMatrixDistributed::fromSharedData(
      this->lhsInverse.getPointer(), processGrid, dim_a, dim_a, parallelConfig.rowBlockSize_,
      parallelConfig.columnBlockSize_);
  this->lhsDistributedInverse.toLocalDataMatrix(this->lhsInverse);
  this->lhsDistributed = DataMatrixDistributed::fromSharedData(
      this->lhsMatrix.data(), processGrid, lhsMatrix.getNrows(), lhsMatrix.getNcols(),
      parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  free(QT);
  free(INV);
  return;
#endif /* USE_SCALAPACK */
}

sgpp::datadriven::MatrixDecompositionType DBMatOfflineOrthoAdapt::getDecompositionType() {
  return sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;
}
}  // namespace datadriven
}  // namespace sgpp
