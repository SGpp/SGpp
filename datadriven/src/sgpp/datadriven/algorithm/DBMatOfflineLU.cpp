// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/base/tools/StringTokenizer.hpp>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

DBMatOfflineLU::DBMatOfflineLU() : DBMatOfflineGE{}, permutation{nullptr} {}

DBMatOfflineLU::DBMatOfflineLU(const DBMatOfflineLU& rhs)
    : DBMatOfflineGE(rhs), permutation(nullptr) {
  size_t gridSize = rhs.lhsMatrix.getNrows();
  permutation = std::unique_ptr<gsl_permutation>{gsl_permutation_alloc(gridSize)};
  gsl_permutation_memcpy(permutation.get(), rhs.permutation.get());
}

DBMatOfflineLU& DBMatOfflineLU::operator=(const DBMatOfflineLU& rhs) {
  size_t gridSize = rhs.lhsMatrix.getNrows();
  DBMatOffline::operator=(rhs);
  permutation = std::unique_ptr<gsl_permutation>{gsl_permutation_alloc(gridSize)};
  gsl_permutation_memcpy(permutation.get(), rhs.permutation.get());

  return *this;
}

DBMatOffline* DBMatOfflineLU::clone() const { return new DBMatOfflineLU{*this}; }

bool DBMatOfflineLU::isRefineable() { return false; }

void DBMatOfflineLU::decomposeMatrix(
    const RegularizationConfiguration& regularizationConfig,
    const DensityEstimationConfiguration& densityEstimationConfig) {
  if (isConstructed) {
    if (isDecomposed) {
      // Already decomposed => Do nothing
      return;
    } else {
      size_t n = lhsMatrix.getNrows();
      gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                                n);  // Create GSL matrix view for decomposition
      permutation =
          std::unique_ptr<gsl_permutation>{gsl_permutation_alloc(n)};  // allocate permutation
      int signum;

      gsl_linalg_LU_decomp(&m.matrix, permutation.get(), &signum);
      isDecomposed = true;
    }

  } else {
    throw base::algorithm_exception("Matrix has to be constructed before it can be decomposed!");
  }
}

DBMatOfflineLU::DBMatOfflineLU(const std::string& fileName)
    : DBMatOfflineGE(), permutation{nullptr} {
  isConstructed = true;
  isDecomposed = true;

  // Read grid size from header (number of rows in lhsMatrix)
  std::ifstream filestream(fileName, std::istream::in);
  // Read configuration
  if (!filestream) {
    throw algorithm_exception("Failed to open File");
  }
  std::string str;
  std::getline(filestream, str);
  filestream.close();

  std::vector<std::string> tokens;
  sgpp::base::StringTokenizer::tokenize(str, ",", tokens);

  auto size = std::stoi(tokens[0]);

  FILE* file = fopen(fileName.c_str(), "rb");
  if (!file) {
    throw algorithm_exception{"Failed to open File"};
  }

  // seek end of first line
  char c = 0;
  while (c != '\n') {
    c = static_cast<char>(fgetc(file));
  }

  // TODO(lettrich) : test if we can do this without copying.
  // Read matrix
  gsl_matrix* matrix;
  matrix = gsl_matrix_alloc(size, size);
  gsl_matrix_fread(file, matrix);

  // read permutation
  permutation = std::unique_ptr<gsl_permutation>{gsl_permutation_alloc(size)};
  gsl_permutation_fread(file, permutation.get());

  fclose(file);

  lhsMatrix = DataMatrix(matrix->data, matrix->size1, matrix->size2);
  gsl_matrix_free(matrix);
}

void DBMatOfflineLU::permuteVector(DataVector& b) {
  if (isDecomposed) {
    gsl_permute(permutation->data, b.getPointer(), 1, b.getSize());
  } else {
    throw algorithm_exception("Matrix was not decomposed yet.");
  }
}

void DBMatOfflineLU::store(const std::string& fileName) {
  // first store header and matrix
  DBMatOffline::store(fileName);
  // then store permutation.
  // c file API needed for GSL
  FILE* outputCFile = fopen(fileName.c_str(), "ab");
  if (!outputCFile) {
    throw algorithm_exception{"cannot open file for writing"};
  }
  gsl_permutation_fwrite(outputCFile, permutation.get());
  fclose(outputCFile);
}

sgpp::datadriven::MatrixDecompositionType DBMatOfflineLU::getDecompositionType() {
  return sgpp::datadriven::MatrixDecompositionType::LU;
}
} /* namespace datadriven */
} /* namespace sgpp */
#endif /* USE_GSL */
