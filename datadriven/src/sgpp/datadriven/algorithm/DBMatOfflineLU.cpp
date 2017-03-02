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

namespace sgpp {
namespace datadriven {} /* namespace datadriven */
} /* namespace sgpp */

sgpp::datadriven::DBMatOfflineLU::DBMatOfflineLU(DBMatDensityConfiguration& oc)
    : DBMatOfflineGE(oc) {}

void sgpp::datadriven::DBMatOfflineLU::decomposeMatrix() {
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

void sgpp::datadriven::DBMatOfflineLU::permuteVector(DataVector& b) {
  if (isDecomposed) {
    gsl_permute(permutation->data, b.getPointer(), 1, b.getSize());
  } else {
    throw base::algorithm_exception("Matrix was not decomposed, yet!");
  }
}
