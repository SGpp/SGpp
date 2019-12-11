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
#include <sgpp/base/grid/type/KinkLinearGrid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitKinkLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModLinear.hpp>

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
#include <set>

namespace sgpp {
namespace datadriven {

using sgpp::base::AdaptivityConfiguration;
using sgpp::base::algorithm_exception;
using sgpp::base::data_exception;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridType;
using sgpp::base::operation_exception;
using sgpp::base::OperationMatrix;
using sgpp::base::RegularGridConfiguration;

DBMatOffline::DBMatOffline()
     : lhsMatrix(), isConstructed(false), isDecomposed(false), lhsInverse(), interactions() {
}

DBMatOffline::DBMatOffline(const DBMatOffline& rhs)
    : lhsMatrix(rhs.lhsMatrix),
      isConstructed(rhs.isConstructed),
      isDecomposed(rhs.isDecomposed),
      lhsInverse(rhs.lhsInverse),
      interactions(rhs.interactions) {}

DBMatOffline& sgpp::datadriven::DBMatOffline::operator=(const DBMatOffline& rhs) {
  if (&rhs == this) {
    return *this;
  }
  lhsMatrix = rhs.lhsMatrix;
  isConstructed = rhs.isConstructed;
  isDecomposed = rhs.isDecomposed;
  lhsInverse = rhs.lhsInverse;
  interactions = rhs.interactions;
  return *this;
}

DBMatOffline::DBMatOffline(const std::string& filepath)
    : lhsMatrix(), isConstructed(true), isDecomposed(true), lhsInverse() {
  // Parse the interactions
  parseInter(filepath, interactions);

  // Parsing of lhsMatrix will be done in subclass implementations
}

DataMatrix& DBMatOffline::getDecomposedMatrix() {
  if (isDecomposed) {
    return lhsMatrix;
  } else {
    throw data_exception("Matrix was not decomposed yet");
  }
}

DataMatrix& DBMatOffline::getInverseMatrix() { return this->lhsInverse; }

DataMatrixDistributed& DBMatOffline::getDecomposedMatrixDistributed() {
#ifdef USE_SCALAPACK
  if (isDecomposed) {
    return lhsDistributed;
  } else {
    throw data_exception("Matrix was not decomposed yet");
  }
#else
  throw sgpp::base::not_implemented_exception("Build without ScaLAPACK");
#endif
}

DataMatrixDistributed& DBMatOffline::getDecomposedInverseDistributed() {
#ifdef USE_SCALAPACK
  if (isDecomposed) {
    return this->lhsDistributedInverse;
  } else {
    throw sgpp::base::algorithm_exception(
        "In DBMatOffline::getDecomposedInverseDistributed:\nlhsMatrix wasn't even decomposed yet.");
  }
#else
  throw sgpp::base::not_implemented_exception("Build without Scalapack");
#endif
}

void DBMatOffline::syncDistributedDecomposition(std::shared_ptr<BlacsProcessGrid> processGrid,
                                                const ParallelConfiguration& parallelConfig) {
#ifdef USE_SCALAPACK
  if (isDecomposed) {
    lhsDistributed = DataMatrixDistributed::fromSharedData(
        lhsMatrix.data(), processGrid, lhsMatrix.getNrows(), lhsMatrix.getNcols(),
        parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);
  } else {
    throw data_exception(
        "In DBMatOffline::syncDistributedDecomposition\nCan't sync, because lhsMatrix "
        "was not decomposed yet");
  }
#endif
  // no action needed without scalapack
}

void DBMatOffline::syncDistributedInverse(std::shared_ptr<BlacsProcessGrid> processGrid,
                                          const ParallelConfiguration& parallelConfig) {
#ifdef USE_SCALAPACK
  if (isDecomposed) {
    lhsDistributedInverse = DataMatrixDistributed::fromSharedData(
        lhsInverse.data(), processGrid, lhsInverse.getNrows(), lhsInverse.getNcols(),
        parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);
  } else {
    throw data_exception(
        "In DBMatOffline::syncDistributedInverse\nCan't sync inverse matrix, because lhsMatrix was "
        "not decomposed yet.");
  }
#endif
  // no action needed without scalapack
}

void DBMatOffline::buildMatrix(Grid* grid,
                               const RegularizationConfiguration& regularizationConfig) {
  if (isConstructed) {  // Already constructed, do nothing
    return;
  }
  // check if grid was created
  if (grid == nullptr) {
    throw algorithm_exception("DBMatOffline: grid was not initialized");
  }

  size_t size = grid->getStorage().getSize();  // Size of the (quadratic) matrices A and C

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
  }

  // Write configuration
  std::ofstream outputFile(fileName, std::ofstream::out);

  if (!outputFile) {
    throw algorithm_exception{"cannot open file for writing"};
  }

  std::string inter = "," + std::to_string(interactions.size());
  for (const std::set<size_t>& i : interactions) {
    inter.append("," + std::to_string(i.size()));
    for (size_t j : i) {
      inter.append("," + std::to_string(j));
    }
  }

  outputFile << lhsMatrix.getNrows() << "," << lhsMatrix.getNcols() << ","
             << static_cast<int>(getDecompositionType()) << inter << "\n";
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
  std::cout << "Stored " << lhsMatrix.getNrows() << "x" << lhsMatrix.getNcols() << " matrix"
            << std::endl;
#else
  throw base::not_implemented_exception("built without GSL");
#endif /* USE_GSL */
}

void DBMatOffline::decomposeMatrixParallel(RegularizationConfiguration& regularizationConfig,
                                           DensityEstimationConfiguration& densityEstimationConfig,
                                           std::shared_ptr<BlacsProcessGrid> processGrid,
                                           const ParallelConfiguration& parallelConfig) {
  throw sgpp::base::algorithm_exception(
      "called ::decomposeMatrixParallel from parent object,\nmeaning decomposition type is not "
      "supported");
}

void DBMatOffline::compute_inverse() {
  throw sgpp::base::algorithm_exception(
      "called ::compute_inverse from parent object,\nmeaning decomposition type is not supported");
}

void DBMatOffline::compute_inverse_parallel(std::shared_ptr<BlacsProcessGrid> processGrid,
                                            const ParallelConfiguration& parallelConfig) {
  throw sgpp::base::algorithm_exception(
      "called ::compute_inverse_parallel from parent object,\nmeaning decomposition type is not "
      "supported");
}

void DBMatOffline::printMatrix() {
  if (isDecomposed) {
    std::cout << "Size: " << lhsMatrix.getNrows() << " , " << lhsMatrix.getNcols() << "\n"
              << lhsMatrix.toString();
  } else {
    throw data_exception("Matrix was not decomposed yet");
  }
}

void DBMatOffline::compute_L2_refine_vectors(DataMatrix* mat_refine, Grid* grid, size_t newPoints) {
  size_t j_start = grid->getSize() - newPoints;

  // note: switch-case when adding support for other gridTypes
  if (grid->getType() == sgpp::base::GridType::Linear) {
    auto opLTwoLin = new sgpp::pde::OperationMatrixLTwoDotExplicitLinear();
    opLTwoLin->buildMatrixWithBounds(mat_refine, grid, 0, 0, j_start, 0);
  } else if (grid->getType() == sgpp::base::GridType::ModLinear) {
    auto opLTwoModLin = new sgpp::pde::OperationMatrixLTwoDotExplicitModLinear();
    opLTwoModLin->buildMatrixWithBounds(mat_refine, grid, 0, 0, j_start, 0);
  } else if (grid->getType() == sgpp::base::GridType::KinkLinear) {
    auto opLTwoKinkLin = new sgpp::pde::OperationMatrixLTwoDotExplicitKinkLinear();
    opLTwoKinkLin->buildMatrixWithBounds(mat_refine, grid, 0, 0, j_start, 0);
  } else {
    throw algorithm_exception(
        "in DBMatOffline::compute_L2_refine_vectors, gridType is not supported.");
  }
}

void sgpp::datadriven::DBMatOffline::parseInter(
    const std::string& fileName, std::set<std::set<size_t>>& interactions) const {
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

  for (size_t i = 4; i < tokens.size(); i += std::stoi(tokens[i]) + 1) {
    std::set<size_t> tmp;
    for (size_t j = 1; j <= std::stoul(tokens[i]); j++) {
      tmp.insert(std::stoi(tokens[i + j]));
    }
    interactions.insert(tmp);
  }

  std::cout << interactions.size() << std::endl;
}

size_t DBMatOffline::getGridSize() { return lhsMatrix.getNrows(); }

sgpp::base::DataMatrix& DBMatOffline::getLhsMatrix_ONLY_FOR_TESTING() { return this->lhsMatrix; }

}  // namespace datadriven
}  // namespace sgpp
