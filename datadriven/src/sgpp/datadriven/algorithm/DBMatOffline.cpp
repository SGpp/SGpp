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

DBMatOffline::DBMatOffline(const DBMatDensityConfiguration& oc)
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
    permutation = std::unique_ptr<gsl_permutation>{gsl_permutation_alloc(size)};
    gsl_permutation_fread(f, permutation.get());
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
    permutation =
        std::unique_ptr<gsl_permutation>{gsl_permutation_alloc(old.grid->getStorage().getSize())};
    gsl_permutation_memcpy(permutation.get(), old.permutation.get());
  }

  // Initialize Grid
  InitializeGrid();
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
    gsl_permutation_fwrite(f, permutation.get());
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

}  // namespace datadriven
}  // namespace sgpp

void sgpp::datadriven::DBMatOffline::buildMatrix() {
  if (isConstructed) {  // Already constructed, do nothing
    return;
  }

  size_t size;

  InitializeGrid();

  // check if grid was created
  if (grid == nullptr) {
    throw application_exception("DBMatOffline: grid was not initialized");
  }

  size = grid->getStorage().getSize();  // Size of the (quadratic) matrices A and C

  // Construct matrix A
  lhsMatrix = DataMatrix(size, size);

  std::unique_ptr<OperationMatrix> op(
      sgpp::op_factory::createOperationLTwoDotExplicit(&lhsMatrix, *grid));
  isConstructed = true;
}

#endif /* USE_GSL */
