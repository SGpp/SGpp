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

DBMatOffline::DBMatOffline(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::base::AdpativityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig)
    : gridConfig(gridConfig), adaptivityConfig(adaptivityConfig),
    regularizationConfig(regularizationConfig), densityEstimationConfig(densityEstimationConfig),
    lhsMatrix(), isConstructed(false), isDecomposed(false), grid(nullptr) {
  interactions = std::vector<std::vector<size_t>>();
}

DBMatOffline::DBMatOffline()
    : gridConfig(), adaptivityConfig(), regularizationConfig(), densityEstimationConfig(),
    lhsMatrix(), isConstructed(false), isDecomposed(false), grid(nullptr) {
  interactions = std::vector<std::vector<size_t>>();
}

DBMatOffline::DBMatOffline(const DBMatOffline& rhs)
    : gridConfig(rhs.gridConfig),
      adaptivityConfig(rhs.adaptivityConfig),
      regularizationConfig(rhs.regularizationConfig),
      densityEstimationConfig(rhs.densityEstimationConfig),
      lhsMatrix(rhs.lhsMatrix),
      isConstructed(rhs.isConstructed),
      isDecomposed(rhs.isDecomposed),
      grid(nullptr),
      interactions(rhs.interactions) {

  if (rhs.grid != nullptr) {
    grid = std::unique_ptr<Grid>{rhs.grid->clone()};
  }
}

DBMatOffline& sgpp::datadriven::DBMatOffline::operator=(const DBMatOffline& rhs) {
  if (&rhs == this) {
    return *this;
  }

  gridConfig = rhs.gridConfig;
  adaptivityConfig = rhs.adaptivityConfig;
  regularizationConfig = rhs.regularizationConfig;
  densityEstimationConfig = rhs.densityEstimationConfig;
  lhsMatrix = rhs.lhsMatrix;
  isConstructed = rhs.isConstructed;
  isDecomposed = rhs.isDecomposed;
  interactions = rhs.interactions;

  if (rhs.grid != nullptr) {
    grid = std::unique_ptr<Grid>{rhs.grid->clone()};
  }
  return *this;
}

DBMatOffline::DBMatOffline(const std::string& fileName)
    : gridConfig(), adaptivityConfig(), regularizationConfig(), densityEstimationConfig(),
    lhsMatrix(), isConstructed(true), isDecomposed(true), grid(nullptr) {
  parseConfig(fileName, gridConfig, adaptivityConfig,
              regularizationConfig, densityEstimationConfig);
  interactions = std::vector<std::vector<size_t>>();
  parseInter(fileName, interactions);
  std::cout << "Setting up Grid..." << std::endl;
  InitializeGrid();
  std::cout << "Grid set up! Start reading Matrix" << std::endl;
}

sgpp::base::GeneralGridConfiguration& DBMatOffline::getGridConfig() {
    return gridConfig;
}
sgpp::base::AdpativityConfiguration& DBMatOffline::getAdaptivityConfig() {
    return adaptivityConfig;
}
sgpp::datadriven::RegularizationConfiguration& DBMatOffline::getRegularizationConfig() {
    return regularizationConfig;
}
sgpp::datadriven::DensityEstimationConfiguration& DBMatOffline::getDensityEstimationConfig() {
    return densityEstimationConfig;
}

DataMatrix& DBMatOffline::getDecomposedMatrix() {
  if (isDecomposed) {
    return lhsMatrix;
  } else {
    throw data_exception("Matrix was not decomposed yet");
  }
}

Grid& DBMatOffline::getGrid() { return *grid; }

void DBMatOffline::InitializeGrid() {
  if (gridConfig.type_ == GridType::ModLinear) {
    grid = std::unique_ptr<Grid>{Grid::createModLinearGrid(gridConfig.dim_)};
  } else if (gridConfig.type_ == GridType::Linear) {
    grid = std::unique_ptr<Grid>{Grid::createLinearGrid(gridConfig.dim_)};
  } else {
    throw algorithm_exception("LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  if (interactions.size() == 0) {
    grid->getGenerator().regular(gridConfig.level_);
  } else {
    grid->getGenerator().regularInter(gridConfig.level_, interactions, 0.0);
  }
  std::cout << "Initialized Grid has " << grid->getSize() << "Gridpoints." << std::endl;
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

  std::string inter = "," + std::to_string(interactions.size());
  for (std::vector<size_t> i : interactions) {
    inter.append("," + std::to_string(i.size()));
    for (size_t j : i) {
      inter.append("," + std::to_string(j));
    }
  }

  outputFile << static_cast<int>(gridConfig.type_) << "," << gridConfig.dim_ << ","
             << gridConfig.level_ << "," << static_cast<int>(regularizationConfig.type_) << ","
             << std::setprecision(12) << regularizationConfig.lambda_ << ","
             << static_cast<int>(densityEstimationConfig.decomposition_) << inter << "\n";
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

void sgpp::datadriven::DBMatOffline::parseConfig(
    const std::string& fileName,
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdpativityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) const {
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

  gridConfig.type_ = static_cast<GridType>(std::stoi(tokens[0]));
  gridConfig.dim_ = std::stoi(tokens[1]);
  gridConfig.level_ = std::stoi(tokens[2]);
  regularizationConfig.type_ = static_cast<RegularizationType>(std::stoi(tokens[3]));
  regularizationConfig.lambda_ = std::stof(tokens[4]);
  densityEstimationConfig.decomposition_ =
  static_cast<MatrixDecompositionType>(std::stoi(tokens[5]));
}

void sgpp::datadriven::DBMatOffline::parseInter(const std::string& fileName,
    std::vector<std::vector<size_t>>& interactions) const {
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

  for (size_t i = 7; i < tokens.size(); i+= std::stoi(tokens[i])+1) {
    std::vector<size_t> tmp = std::vector<size_t>();
    for (size_t j = 1; j <= std::stoul(tokens[i]); j++) {
      tmp.push_back(std::stoi(tokens[i+j]));
    }
    interactions.push_back(tmp);
  }

  std::cout << interactions.size() << std::endl;
}

void sgpp::datadriven::DBMatOffline::setInter(std::vector<std::vector <size_t>> inter) {
  interactions = inter;
}
sgpp::base::DataMatrix& DBMatOffline::getLhsMatrix_ONLY_FOR_TESTING() { return this->lhsMatrix; }

}  // namespace datadriven
}  // namespace sgpp
