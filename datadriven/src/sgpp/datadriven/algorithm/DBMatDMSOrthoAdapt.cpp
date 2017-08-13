// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>

namespace sgpp {
namespace datadriven {

void DBMatDMSOrthoAdapt::solve(sgpp::base::DataMatrix& T_inv, sgpp::base::DataMatrix& Q,
                               sgpp::base::DataMatrix& B, sgpp::base::DataVector& b,
                               sgpp::base::DataVector& alpha) {
  // assert dimensions
  bool dimension_check = Q.getNrows() == b.getSize() && B.getNcols() == alpha.getSize();
  if (!dimension_check) {
    throw sgpp::base::data_exception("Dimensions of b or alpha do not match the system");
  }

  // stores interim values
  sgpp::base::DataVector buffer(alpha.getSize(), 0.0);

  // start with Q^t * x, doing this manually until sgpp supports transpose operations
  for (size_t i = 0; i < Q.getNcols(); i++) {
    double acc = 0.0;
    for (size_t j = 0; j < Q.getNrows(); j++) {
      acc += Q.get(j, i) * b.get(j);  // note: i and j are switched in Q
    }
    alpha.set(i, acc);
  }

  // T_inv * Q^t * x
  T_inv.mult(alpha, buffer);

  // Q * T_inv * Q^t * x
  Q.mult(buffer, alpha);

  // B * b
  B.mult(b, buffer);

  // Q*T_inv*Q^t*b + B*b
  alpha.add(buffer);
}
}  // namespace datadriven
}  // namespace sgpp
