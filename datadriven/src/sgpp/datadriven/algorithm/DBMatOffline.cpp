// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#ifdef USE_GSL
#include <gsl/gsl_matrix_double.h>
#endif /* USE_GSL */

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
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
using sgpp::base::algorithm_exception;
using sgpp::base::data_exception;
using sgpp::base::OperationMatrix;

DBMatOffline::DBMatOffline(const DBMatDensityConfiguration& oc)
    : config(oc), lhsMatrix(), isConstructed(false), isDecomposed(false), grid(nullptr) {}

DBMatOffline::DBMatOffline()
    : config(), lhsMatrix(), isConstructed(false), isDecomposed(false), grid(nullptr) {}

DBMatOffline::DBMatOffline(const DBMatOffline& rhs)
    : config(rhs.config),
      lhsMatrix(rhs.lhsMatrix),
      isConstructed(rhs.isConstructed),
      isDecomposed(rhs.isDecomposed),
      grid(nullptr) {
  if (rhs.grid != nullptr) {
    grid = std::unique_ptr<Grid>{rhs.grid->clone()};
  }
}

DBMatOffline& sgpp::datadriven::DBMatOffline::operator=(const DBMatOffline& rhs) {
  if (&rhs == this) {
    return *this;
  }

  config = rhs.config;
  lhsMatrix = rhs.lhsMatrix;
  isConstructed = rhs.isConstructed;
  isDecomposed = rhs.isDecomposed;

  if (rhs.grid != nullptr) {
    grid = std::unique_ptr<Grid>{rhs.grid->clone()};
  }
  return *this;
}

DBMatOffline::DBMatOffline(const std::string& fileName)
    : config(), lhsMatrix(), isConstructed(true), isDecomposed(true), grid(nullptr) {
  parseConfig(fileName, config);
  InitializeGrid();
}

DBMatDensityConfiguration& DBMatOffline::getConfig() { return config; }

DataMatrix& DBMatOffline::getDecomposedMatrix() {
  if (isDecomposed) {
    return lhsMatrix;
  } else {
    throw data_exception("Matrix was not decomposed yet");
  }
}

Grid& DBMatOffline::getGrid() { return *grid; }

void DBMatOffline::InitializeGrid() {
  if (config.grid_type_ == GridType::ModLinear) {
    grid = std::unique_ptr<Grid>{Grid::createModLinearGrid(config.grid_dim_)};
  } else if (config.grid_type_ == GridType::Linear) {
    grid = std::unique_ptr<Grid>{Grid::createLinearGrid(config.grid_dim_)};
  } else {
    throw algorithm_exception("LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
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
  if (grid == nullptr) {
    throw algorithm_exception("DBMatOffline: grid was not initialized");
  }

  size = grid->getStorage().getSize();  // Size of the (quadratic) matrices A and C

  // Construct matrix A
  lhsMatrix = DataMatrix(size, size);

  std::unique_ptr<OperationMatrix> op(
      op_factory::createOperationLTwoDotExplicit(&lhsMatrix, *grid));
  isConstructed = true;
}

void DBMatOffline::store(const std::string& fileName) {
#ifdef USE_GSL
  if (!isDecomposed) {
    throw algorithm_exception("Matrix not decomposed yet");
    return;
  }

  // Write configuration
  std::ofstream outputFile(fileName, std::ofstream::out);

  if (!outputFile) {
    throw algorithm_exception{"cannot open file for writing"};
  }

  outputFile << static_cast<int>(config.grid_type_) << "," << config.grid_dim_ << ","
             << config.grid_level_ << "," << static_cast<int>(config.regularization_) << ","
             << std::setprecision(12) << config.lambda_ << ","
             << static_cast<int>(config.decomp_type_) << "\n";
  outputFile.close();

  // write matrix
  // switch to c FILE API for GSL
  FILE* outputCFile = fopen(fileName.c_str(), "ab");
  if (!outputCFile) {
    throw algorithm_exception{"cannot open file for writing"};
  }
  gsl_matrix_view matrixView =
      gsl_matrix_view_array(lhsMatrix.getPointer(), lhsMatrix.getNrows(), lhsMatrix.getNcols());
  gsl_matrix_fwrite(outputCFile, &matrixView.matrix);

  fclose(outputCFile);
#else
  throw base::not_implemented_exception("built withot GSL");
#endif /* USE_GSL */
}

void DBMatOffline::printMatrix() {
  if (isDecomposed) {
    std::cout << "Size: " << lhsMatrix.getNrows() << " , " << lhsMatrix.getNcols() << "\n"
              << lhsMatrix.toString();
  } else {
    throw data_exception("Matrix was not decomposed yet");
  }
}

void sgpp::datadriven::DBMatOffline::parseConfig(const std::string& fileName,
                                                 DBMatDensityConfiguration& config) const {
  std::ifstream file(fileName, std::istream::in);
  // Read configuration
  if (!file) {
    throw algorithm_exception("Failed to open File");
  }
  std::string str;
  std::getline(file, str);
  file.close();

  std::vector<std::string> tokens;
  StringTokenizer::tokenize(str, ",", tokens);

  config.grid_type_ = static_cast<GridType>(std::stoi(tokens[0]));
  config.grid_dim_ = std::stoi(tokens[1]);
  config.grid_level_ = std::stoi(tokens[2]);
  config.regularization_ = static_cast<RegularizationType>(std::stoi(tokens[3]));
  config.lambda_ = std::stof(tokens[4]);
  config.decomp_type_ = static_cast<DBMatDecompostionType>(std::stoi(tokens[5]));
}

}  // namespace datadriven
}  // namespace sgpp
