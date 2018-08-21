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

DBMatOffline::DBMatOffline() : lhsMatrix(), isConstructed(false), isDecomposed(false) {
  interactions = std::vector<std::vector<size_t>>();
}

DBMatOffline::DBMatOffline(const DBMatOffline& rhs)
    : lhsMatrix(rhs.lhsMatrix),
      isConstructed(rhs.isConstructed),
      isDecomposed(rhs.isDecomposed),
      interactions(rhs.interactions) { }

DBMatOffline& sgpp::datadriven::DBMatOffline::operator=(const DBMatOffline& rhs) {
  if (&rhs == this) {
    return *this;
  }
  lhsMatrix = rhs.lhsMatrix;
  isConstructed = rhs.isConstructed;
  isDecomposed = rhs.isDecomposed;
  interactions = rhs.interactions;
  return *this;
}

DBMatOffline::DBMatOffline(const std::string& filepath)
    : lhsMatrix(), isConstructed(true), isDecomposed(true) {
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


void DBMatOffline::buildMatrix(Grid* grid, RegularizationConfiguration& regularizationConfig) {
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
  std::cout << "Stored " << lhsMatrix.getNrows() << "x" << lhsMatrix.getNcols() << " matrix" <<
      std::endl;
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

  for (size_t i = 4; i < tokens.size(); i+= std::stoi(tokens[i])+1) {
    std::vector<size_t> tmp = std::vector<size_t>();
    for (size_t j = 1; j <= std::stoul(tokens[i]); j++) {
      tmp.push_back(std::stoi(tokens[i+j]));
    }
    interactions.push_back(tmp);
  }

  std::cout << interactions.size() << std::endl;
}


size_t DBMatOffline::getGridSize() { return lhsMatrix.getNrows(); }

sgpp::base::DataMatrix& DBMatOffline::getLhsMatrix_ONLY_FOR_TESTING() { return this->lhsMatrix; }

}  // namespace datadriven
}  // namespace sgpp
