/*
 * DBMatOfflineLU.cpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>

#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOfflineLU::DBMatOfflineLU(const DBMatDensityConfiguration& oc)
    : DBMatOfflineGE{oc}, permutation{nullptr} {}



    void
    DBMatOfflineLU::decomposeMatrix() {
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

void DBMatOfflineLU::permuteVector(DataVector& b) {
  if (isDecomposed) {
    gsl_permute(permutation->data, b.getPointer(), 1, b.getSize());
  } else {
    throw algorithm_exception("Matrix was not decomposed, yet!");
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

} /* namespace datadriven */
} /* namespace sgpp */
