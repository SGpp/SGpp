// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;

sgpp::datadriven::DBMatOnlineDELU::DBMatOnlineDELU(DBMatOffline& offline, Grid& grid, double lambda,
    double beta)
    : DBMatOnlineDE{offline, grid, lambda, beta} {}

void sgpp::datadriven::DBMatOnlineDELU::solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool do_cv) {
  DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();

  // Solve the system:
  alpha = DataVector(lhsMatrix.getNcols());
  DBMatDMSBackSub lusolver;
  lusolver.solve(lhsMatrix, alpha, b);
}

} /* namespace datadriven */
} /* namespace sgpp */
#endif /*USE_GSL*/
