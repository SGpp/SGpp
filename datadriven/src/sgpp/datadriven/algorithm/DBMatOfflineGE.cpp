// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/tools/StringTokenizer.hpp>

#ifdef USE_GSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#endif /* USE_GSL */

#include <list>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;
using sgpp::base::application_exception;
using sgpp::base::DataMatrix;
using sgpp::base::operation_exception;

DBMatOfflineGE::DBMatOfflineGE() : DBMatOffline() {}

sgpp::datadriven::DBMatOfflineGE::DBMatOfflineGE(const std::string& fileName)
    : DBMatOffline{fileName} {
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

#ifdef USE_GSL
  auto size = std::stoi(tokens[0]);

  FILE* file = fopen(fileName.c_str(), "rb");
  if (!file) {
    throw application_exception{"Failed to open File"};
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
  fclose(file);

  lhsMatrix = DataMatrix(matrix->data, matrix->size1, matrix->size2);
  gsl_matrix_free(matrix);
#else
  throw base::not_implemented_exception("built withot GSL");
#endif /* USE_GSL */
}

void DBMatOfflineGE::buildMatrix(Grid* grid,
                                 const RegularizationConfiguration& regularizationConfig) {
  // build matrix
  DBMatOffline::buildMatrix(grid, regularizationConfig);

  // then add regularization term
  auto size = grid->getStorage().getSize();

  // Construct matrix lambda * C (just use identity for C)
  DataMatrix lambdaC(size, size);
  if (regularizationConfig.type_ == RegularizationType::Identity) {
    lambdaC.setAll(0.);
    for (size_t i = 0; i < size; i++) {
      lambdaC.set(i, i, regularizationConfig.lambda_);
    }
  } else {
    throw operation_exception("Unsupported regularization type");
  }

  // Compute A + lambda * C:
  lhsMatrix.add(lambdaC);

  isConstructed = true;
}

} /* namespace datadriven */
} /* namespace sgpp */
