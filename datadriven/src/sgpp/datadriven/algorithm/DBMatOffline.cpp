// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <list>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::Grid;
using sgpp::base::GridType;
using sgpp::base::RegularGridConfiguration;
using sgpp::base::AdpativityConfiguration;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::operation_exception;
using sgpp::base::application_exception;
using sgpp::base::data_exception;
using sgpp::base::OperationMatrix;

DBMatOffline::DBMatOffline(DBMatDensityConfiguration& oc)
    : config(oc),
      lhsMatrix(),
      isConstructed(false),
      isDecomposed(false),
      permutation(nullptr),
      grid(nullptr) {}

DBMatOffline::DBMatOffline()
    : config(),
      lhsMatrix(),
      isConstructed(false),
      isDecomposed(false),
      permutation(nullptr),
      grid(nullptr) {}

DBMatOffline::DBMatOffline(const std::string& fname)
    : config(),
      lhsMatrix(),
      isConstructed(false),
      isDecomposed(false),
      permutation(nullptr),
      grid(nullptr) {
  // std::cout << "START READING MATRIX" << std::endl;
  // SGppStopwatch* myStopwatch = new SGppStopwatch();
  // myStopwatch->start();

  // Read configuration
  FILE* f = fopen(fname.c_str(), "rb");
  std::string line("");
  if (f == nullptr) {
    std::cout << "DBMatOffline: Error opening file " << line << std::endl;
    exit(-1);
  }
  char c = static_cast<char>(fgetc(f));
  line += c;
  while (c != '\n') {
    c = static_cast<char>(fgetc(f));
    line += c;
  }

  std::vector<std::string> tokens;
  std::string delim(",");
  Tokenize(line, tokens, delim);

  // ToDo: full grid not supported
  bool fullgrid = atoi(tokens[0].c_str());
  if ((fullgrid && tokens.size() != 8) || (!fullgrid && tokens.size() != 7)) {
    std::cout << "DBMatOffline: Wrong file format: " << fname.c_str() << std::endl;
    exit(-1);
  }

  GridType grid_type = (GridType)atoi(tokens[1].c_str());

  size_t grid_dim = atoi(tokens[2].c_str());
  int grid_level = atoi(tokens[3].c_str());
  RegularizationType reg = (RegularizationType)atoi(tokens[4].c_str());
  double lambda = atof(tokens[5].c_str());
  DBMatDecompostionType decomp = (DBMatDecompostionType)atoi(tokens[6].c_str());

  RegularGridConfiguration gconf;
  gconf.dim_ = grid_dim;
  gconf.level_ = grid_level;
  gconf.type_ = grid_type;

  base::AdpativityConfiguration adaptivityConfig;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.threshold_ = 0.0;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.percent_ = 0.0;

  config = DBMatDensityConfiguration(gconf, adaptivityConfig, reg, lambda, decomp);

  size_t size;

  // Build grid
  InitializeGrid();

  // check if grid was created
  if (grid == nullptr) return;

  size = grid->getStorage().getSize();  // Size of the (quadratic) matrices A and C

  // Read matrix
  gsl_matrix* a;
  if (decomp == DBMatDecompostionType::DBMatDecompLU) {
    a = gsl_matrix_alloc(size, size);
  } else if (decomp == DBMatDecompostionType::DBMatDecompEigen) {
    a = gsl_matrix_alloc(size + 1, size);
  } else if (decomp == DBMatDecompostionType::DBMatDecompChol) {
    a = gsl_matrix_alloc(size, size);
  } else {
    throw application_exception("Unsupported decomposition type!");
  }

  gsl_matrix_fread(f, a);
  lhsMatrix = DataMatrix(a->data, a->size1, a->size2);
  // Read permutation
  if (decomp == DBMatDecompostionType::DBMatDecompLU) {
    permutation = gsl_permutation_alloc(size);
    gsl_permutation_fread(f, permutation);
  }

  fclose(f);
  gsl_matrix_free(a);

  isConstructed = true;
  isDecomposed = true;

  // std::cout << "Time: " << myStopwatch->stop() << std::endl;
}

DBMatOffline::DBMatOffline(const DBMatOffline& old)
    : config(),
      lhsMatrix(),
      isConstructed(true),
      isDecomposed(true),
      permutation(nullptr),
      grid(nullptr) {
  if (!old.isDecomposed) {
    throw application_exception("Matrix not yet decomposed");
    return;
  }

  // Copy configuration
  RegularGridConfiguration gconf;
  gconf.dim_ = old.config.grid_dim_;
  gconf.level_ = old.config.grid_level_;
  gconf.type_ = old.config.grid_type_;

  AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = old.config.numRefinements_;
  adaptConfig.threshold_ = old.config.ref_threshold_;
  adaptConfig.noPoints_ = old.config.ref_noPoints_;

  config = DBMatDensityConfiguration{gconf, adaptConfig, old.config.regularization_,
                                     old.config.lambda_, old.config.decomp_type_};

  // Copy Matrix
  lhsMatrix = DataMatrix(old.lhsMatrix);

  // Copy Permutation (if existing)
  if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompLU) {
    permutation = gsl_permutation_alloc(old.grid->getStorage().getSize());
    gsl_permutation_memcpy(permutation, old.permutation);
  }

  // Initialize Grid
  InitializeGrid();
}

DBMatOffline::~DBMatOffline() {
  if (permutation != nullptr) delete permutation;
}

DBMatDensityConfiguration& DBMatOffline::getConfig() { return config; }

DataMatrix& DBMatOffline::getDecomposedMatrix() {
  if (isDecomposed) {
    return lhsMatrix;
  } else {
    throw data_exception("Matrix was not decomposed, yet!");
  }
}

Grid& DBMatOffline::getGrid() { return *grid; }

void DBMatOffline::permuteVector(DataVector& b) {
  if (isDecomposed) {
    if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompLU) {
      gsl_permute(permutation->data, b.getPointer(), 1, b.getSize());
    } else {
      // No permutation needed for other decomposition types!
    }
  } else {
    throw data_exception("Matrix was not decomposed, yet!");
  }
}

void DBMatOffline::InitializeGrid() {
  if (config.grid_type_ == GridType::ModLinear) {
    grid = std::unique_ptr<Grid>{Grid::createModLinearGrid(config.grid_dim_)};
  } else if (config.grid_type_ == GridType::Linear) {
    grid = std::unique_ptr<Grid>{Grid::createLinearGrid(config.grid_dim_)};
  } else {
    throw application_exception(
        "LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  grid->getGenerator().regular(config.grid_level_);
}

void DBMatOffline::buildMatrix() {
  if (isConstructed) {  // Already constructed, do nothing
    return;
  }

  size_t size;

  InitializeGrid();

  // check if grid was created
  if (grid == nullptr) return;

  size = grid->getStorage().getSize();  // Size of the (quadratic) matrices A and C

  // Construct matrix A
  lhsMatrix = DataMatrix(size, size);

  std::unique_ptr<OperationMatrix> op(
      sgpp::op_factory::createOperationLTwoDotExplicit(&lhsMatrix, *grid));

  if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompLU ||
      config.decomp_type_ == DBMatDecompostionType::DBMatDecompChol) {
    // Add regularization term:
    // Construct matrix lambda * C (just use identity for C)
    DataMatrix lambdaC(size, size);
    if (config.regularization_ == RegularizationType::Identity) {
      lambdaC.setAll(0.);
      for (size_t i = 0; i < size; i++) {
        lambdaC.set(i, i, config.lambda_);
      }
    } else {
      throw operation_exception("Unsupported regularization type");
    }

    // Compute A + lambda * C:
    lhsMatrix.add(lambdaC);
  }

  isConstructed = true;
}

void DBMatOffline::decomposeMatrix() {
  if (isConstructed) {
    if (isDecomposed) {
      // Already decomposed => Do nothing
      return;
    } else {
      size_t n = lhsMatrix.getNrows();
      if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompLU) {
        gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                                  n);  // Create GSL matrix view for decomposition
        int* sig = reinterpret_cast<int*>(malloc(sizeof(int)));
        permutation = gsl_permutation_alloc(n);  // allocate permutation

        gsl_linalg_LU_decomp(&m.matrix, permutation, sig);

        delete sig;

      } else if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompEigen) {
        gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                                  n);  // Create GSL matrix view for decomposition

        gsl_matrix* q = gsl_matrix_alloc(n, n);  // Stores the eigenvectors
        gsl_vector* e = gsl_vector_alloc(n);     // Stores the eigenvalues

        gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc(n);
        gsl_eigen_symmv(&m.matrix, e, q, ws);
        gsl_eigen_symmv_free(ws);

        // Create an (n+1)*n matrix to store eigenvalues and -vectors:
        lhsMatrix = DataMatrix(n + 1, n);

        for (size_t r = 0; r < n; r++) {
          for (size_t c = 0; c < n; c++) {
            lhsMatrix.set(r, c, gsl_matrix_get(q, r, c));
          }
        }
        for (size_t c = 0; c < n; c++) {
          lhsMatrix.set(n, c, gsl_vector_get(e, c));
        }

        delete e;
        delete q;
      } else if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompChol) {
        gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                                  n);  // Create GSL matrix view for decomposition
        // Perform Cholesky decomposition
        gsl_linalg_cholesky_decomp(&m.matrix);

        // Isolate lower triangular matrix
        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < n; j++) {
            if (i < j) {
              lhsMatrix.set(i, j, 0);
            }
          }
        }
      } else {
        // No other decomposition types, yet
        return;
      }
      isDecomposed = true;
    }
  } else {
    throw data_exception("Matrix has to be constructed before it can be decomposed!");
  }
}

void DBMatOffline::store(const std::string& fileName) {
  if (!isDecomposed) {
    throw application_exception("Matrix not yet decomposed");
    return;
  }

  // Write configuration
  FILE* f = fopen(fileName.c_str(), "w");
  if (!(f != 0)) {
    std::cout << "libtool: DBMatOffline: Error opening file " << fileName.c_str() << std::endl;
    exit(-1);
  }

  /*fprintf(f, "%d,%d,%d,%d,%d,%.10e,%d", config_->grid_type_, config_->grid_dim_,
          config_->grid_level_, config_->regularization_, config_->lambda_,
          config_->decomp_type_);
  fprintf(f, "\n");
  fclose(f);*/
  // Write Matrix
  f = fopen(fileName.c_str(), "ab");
  gsl_matrix_view m =
      gsl_matrix_view_array(lhsMatrix.getPointer(), lhsMatrix.getNrows(), lhsMatrix.getNcols());
  gsl_matrix_fwrite(f, &m.matrix);

  // Write Permutation (if existing)
  if (config.decomp_type_ == DBMatDecompostionType::DBMatDecompLU) {
    gsl_permutation_fwrite(f, permutation);
  }

  fclose(f);
}

void DBMatOffline::printMatrix() {
  std::cout << "Size: " << lhsMatrix.getNrows() << " , " << lhsMatrix.getNcols() << std::endl;
  for (size_t r = 0; r < lhsMatrix.getNrows(); r++) {
    for (size_t c = 0; c < lhsMatrix.getNcols(); c++) {
      std::cout << lhsMatrix.get(r, c) << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void DBMatOffline::Tokenize(std::string& str, std::vector<std::string>& tokens,
                            std::string& delimiters) {
  /*if (!strcmp(delimiters.c_str(), "")) {*/
  if (!delimiters.compare("")) {
    delimiters = " ";
  }
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void DBMatOffline::choleskyModification(size_t newPoints, std::list<size_t> deletedPoints,
                                        double lambda) {
  // Start coarsening
  // If list 'deletedPoints' is not empty, grid points got removed
  if (deletedPoints.size() > 0) {
    size_t new_size = grid->getSize();
    size_t old_size = new_size - newPoints + deletedPoints.size();

    // 'c' is the threshold to decide what kind of Cholesky modification should
    // be performed to take
    // account for the removed grid points
    double c = static_cast<double>(old_size) / 10;
    // Indicates which permutation to perform
    bool start1 = false;
    size_t index_coarse = 0;

    // If index of removed grid point <= c-> Permutation towards first
    // row/column -> job = 1
    size_t coarseCount_1 = 0;
    // If index of removed grid point > c -> Permutation towards last row/column
    // -> job = 2
    size_t coarseCount_2 = 0;

    for (std::list<size_t>::reverse_iterator it = deletedPoints.rbegin();
         it != deletedPoints.rend(); it++) {
      // Indicates current row/column to permute
      index_coarse = *it + 1 + coarseCount_1;

      if (index_coarse > c && start1 == false) {
        // Is accessed if index is larger than 'c'
        choleskyPermutation(index_coarse, old_size - coarseCount_2, 2);
        coarseCount_2++;
      } else if (start1 == false) {
        // Is accessed the first 'index_coarse' is less than 'c'
        start1 = true;
        // Delete last 'coarseCount_2' rows/columns of current Cholesky factor
        lhsMatrix.resizeQuadratic(old_size - coarseCount_2);
      }

      if (start1 == true) {
        // Is accessed if 'index_coarse' is less than 'c'
        choleskyPermutation(1 + coarseCount_1, index_coarse, 1);
        coarseCount_1++;
      }
    }

    if (coarseCount_1 > 0) {
      // If some indices have been less than 'c'
      DataMatrix* update_matrix = new DataMatrix(lhsMatrix.getNrows(), lhsMatrix.getNcols());
      update_matrix->copyFrom(lhsMatrix);
      // Resize copy of current Cholesky factor to
      // receive a matrix of rank one update vectors
      update_matrix->resizeToSubMatrix(coarseCount_1 + 1, 1, lhsMatrix.getNrows(), coarseCount_1);
      // Resize current Cholesky factor to required submatrix
      // for necessary rank one updates
      lhsMatrix.resizeToSubMatrix(coarseCount_1 + 1, coarseCount_1 + 1, lhsMatrix.getNrows(),
                                  lhsMatrix.getNrows());
      DataVector* temp_col = new DataVector(update_matrix->getNrows());

      // 'coarseCount_1' many rank one updates based on the columns of
      // 'update_matrix' are performed
      DBMatDMSChol cholsolver;
      for (size_t i = 0; i < coarseCount_1; i++) {
        update_matrix->getColumn(i, *temp_col);
        cholsolver.choleskyUpdate(lhsMatrix, temp_col, false);
      }
      delete update_matrix;
      delete temp_col;
    } else {
      // If no indices have been less than 'c'
      lhsMatrix.resizeQuadratic(old_size - coarseCount_2);
    }
  }

  // Start refinement
  if (newPoints > 0) {
    size_t gridSize = grid->getStorage().getSize();
    size_t gridDim = grid->getStorage().getDimension();

    // DataMatrix to collect vectors to append
    DataMatrix* mat_refine = new DataMatrix(gridSize, newPoints);

    DataMatrix level(gridSize, gridDim);
    DataMatrix index(gridSize, gridDim);

    grid->getStorage().getLevelIndexArraysForEval(level, index);
    double lambda_conf = lambda;
    // Loop to calculate all L2-products of added points based on the
    // hat-function as basis function
    for (size_t i = 0; i < gridSize; i++) {
      for (size_t j = gridSize - newPoints; j < gridSize; j++) {
        double res = 1;
        for (size_t k = 0; k < gridDim; k++) {
          double lik = level.get(i, k);
          double ljk = level.get(j, k);
          double iik = index.get(i, k);
          double ijk = index.get(j, k);

          if (lik == ljk) {
            if (iik == ijk) {
              // Use formula for identical ansatz functions:
              res *= 2 / lik / 3;
            } else {
              // Different index, but same level => ansatz functions do not
              // overlap:
              res = 0.;
              break;
            }
          } else {
            if (std::max((iik - 1) / lik, (ijk - 1) / ljk) >=
                std::min((iik + 1) / lik, (ijk + 1) / ljk)) {
              // Ansatz functions do not not overlap:
              res = 0.;
              break;
            } else {
              // Use formula for different overlapping ansatz functions:
              if (lik > ljk) {  // Phi_i_k is the "smaller" ansatz function
                double diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
                double temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
                temp_res *= ljk;
                temp_res = (1 - temp_res) / lik;
                res *= temp_res;
              } else {  // Phi_j_k is the "smaller" ansatz function
                double diff = (ijk / ljk) - (iik / lik);  // x_j_k - x_i_k
                double temp_res = fabs(diff - (1 / ljk)) + fabs(diff + (1 / ljk)) - fabs(diff);
                temp_res *= lik;
                temp_res = (1 - temp_res) / ljk;
                res *= temp_res;
              }
            }
          }
        }
        // The new Rows/Cols are stored in mat_refine

        // add current lambda to lower diagonal elements of mat_refine
        if (i == j) {
          mat_refine->set(i, j - gridSize + newPoints, res + lambda_conf);
        } else {
          mat_refine->set(i, j - gridSize + newPoints, res);
        }
      }
    }

    // Resize Cholesky factor to new 'gridSize' before 'choleskyAddPoint' is
    // applied
    // in order to save runtime
    this->lhsMatrix.resizeQuadratic(gridSize);

    // Now its time to call 'choleskyAddPoint''countNewGridPoints' often
    DataVector* temp_col = new DataVector(gridSize);
    for (size_t j = gridSize - newPoints; j < gridSize; j++) {
      temp_col->resizeZero(gridSize);
      mat_refine->getColumn(j - gridSize + newPoints, *temp_col);
      temp_col->resizeZero(j + 1);
      choleskyAddPoint(temp_col, j);
    }

    delete temp_col;
  }

  return;
}

void DBMatOffline::choleskyAddPoint(DataVector* newCol, size_t size) {
  if (!isDecomposed) {
    throw data_exception("Matrix was not decomposed, yet!");
  }

  DataMatrix& mat = lhsMatrix;
  // Size of provided memory for Cholesky factor,
  // because the allocations take place in 'choleskyModifications'
  size_t size_full = mat.getNrows();
  // Size of Cholesky factor after adding 'newCol'
  size_t size_up = newCol->getSize();

  if (size_up != (size + 1)) {
    throw data_exception(
        "Size of update vector newCol needs to be 1 dim larger than the "
        "underlying decomposed matrix!");
  }

  // Create GSL matrix view for update procedures
  gsl_matrix_view m_full = gsl_matrix_view_array(mat.getPointer(), size_full, size_full);
  // Access submatrx since mat_ has already expanded to save alocations
  // procedures
  gsl_matrix_view m = gsl_matrix_submatrix(&m_full.matrix, 0, 0, size, size);

  gsl_vector_view vvec = gsl_vector_view_array(newCol->getPointer(), size_up);
  gsl_vector* wkvec_a = gsl_vector_calloc(size);

  // Extract newCol(0:size-1) = a
  gsl_vector_view c = gsl_vector_subvector(&vvec.vector, 0, size);

  // Solve system a = Lc
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &m.matrix, &c.vector);
  gsl_blas_dcopy(&c.vector, wkvec_a);

  double phi;
  gsl_blas_ddot(wkvec_a, &c.vector, &phi);

  // Compute d = sqrt(newCol[last] - phi)
  double last = gsl_vector_get(&vvec.vector, size_up - 1);
  phi = last - phi;
  // Ensure 'phi' is larger 0
  if (phi <= 0) {
    throw data_exception("Resulting matrix is at least not numerical positive definite!");
  } else {
    phi = sqrt(phi);
  }

  // Modify Choleskyfactor m (L)

  // DataVector full of zeros
  DataVector* zeros = new DataVector(size_full, 0.0);

  // Add 'newCol' and 'zeros' to current Cholesky factor 'lhsMatrix_'
  mat.setColumn(size, *zeros);
  newCol->set(size_up - 1, phi);
  newCol->resizeZero(size_full);
  mat.setRow(size, *newCol);

  gsl_vector_free(wkvec_a);
  delete zeros;

  return;
}

void DBMatOffline::choleskyPermutation(size_t k, size_t l, size_t job) {
  if (!isDecomposed) {
    throw data_exception("Matrix was not decomposed, yet!");
  }

  DataMatrix& mat = lhsMatrix;
  size_t size = mat.getNrows();

  if (k > l) {
    throw data_exception("l needs to be larger than k");
  } else if (l > size) {
    throw data_exception("l needs to be smaller than the column size of the matrix");
  } else if (l == k) {
    return;
  }

  // Determine amount of necessary Givens rotations and rows to swap
  size_t count_perm = l - k;
  // Create GSL matrix view for update procedures
  gsl_matrix_view m = gsl_matrix_view_array(mat.getPointer(), size, size);

  // Define and declare Workingvector, Cosine- and Sinevector
  gsl_vector* svec = gsl_vector_calloc(count_perm);
  gsl_vector* cvec = gsl_vector_calloc(count_perm);

  double* givens_zero;
  double* tbuff = m.matrix.data;

  if (job == 2) {
    // Permute upper triangular - job = 2        => left circular shift
    // 1,...,k-1,k,k+1, ..., l-1,l,l+1, ..,size  => 1,...,k-1,k+1, ...,
    // l-1,l,k,l+1,..., size
    for (size_t i = k - 1; i < l - 1; i++) {
      gsl_matrix_swap_rows(&m.matrix, i, i + 1);
    }

    tbuff += (k - 1) * size + (k - 1);
    for (size_t j = 1; j <= count_perm; j++) {
      givens_zero = gsl_matrix_ptr(&m.matrix, k + j - 2, k + j - 1);
      // Givensrotation in (k+i-1,k+i)-Plane
      gsl_blas_drotg(tbuff, givens_zero, cvec->data + j - 1, svec->data + j - 1);
      // Access columns to modify via Givens rotation
      gsl_vector_view diag_sub =
          gsl_matrix_subcolumn(&m.matrix, k + j - 2, k + j - 1, size - k - j + 1);
      gsl_vector_view low_diag_sub =
          gsl_matrix_subcolumn(&m.matrix, k + j - 1, k + j - 1, size - k - j + 1);
      gsl_matrix_set(&m.matrix, k + j - 2, k + j - 1, 0.0);
      // Apply Givens rotation
      gsl_blas_drot(&diag_sub.vector, &low_diag_sub.vector, cvec->data[j - 1], svec->data[j - 1]);
      tbuff += (size + 1);
    }
  } else if (job == 1) {
    // Permute upper triangular - job = 1       => right circular shift
    // 1,...,k-1,k,k+1, ..., l-1,l,l+1,...size  => 1,...,k-1,l,k,k+1, ...,
    // l-1,l+1,...size
    for (size_t i = l - 1; i > k - 1; i--) {
      gsl_matrix_swap_rows(&m.matrix, i, i - 1);
    }

    tbuff += (k - 1) * size + (l - 2);
    for (size_t j = 1; j <= count_perm; j++) {
      givens_zero = gsl_matrix_ptr(&m.matrix, k - 1, l - j);
      // Givensrotation in (l-i,l-i+1)-Plane
      gsl_blas_drotg(tbuff, givens_zero, cvec->data + j - 1, svec->data + j - 1);
      // Access columns to modify via Givens rotation
      gsl_vector_view diag_sub = gsl_matrix_subcolumn(&m.matrix, l - j - 1, l - j, size - l + j);
      gsl_vector_view low_diag_sub = gsl_matrix_subcolumn(&m.matrix, l - j, l - j, size - l + j);
      gsl_matrix_set(&m.matrix, k - 1, l - j, 0.0);
      // Apply Givens rotation
      gsl_blas_drot(&diag_sub.vector, &low_diag_sub.vector, cvec->data[j - 1], svec->data[j - 1]);
      tbuff -= 1;
    }
  }

  gsl_vector_free(svec);
  gsl_vector_free(cvec);

  return;
}

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
