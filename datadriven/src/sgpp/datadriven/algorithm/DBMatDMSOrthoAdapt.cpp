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
  if (B.getNcols() != 1 && B.getNcols() != b.getSize()) {
    throw sgpp::base::data_exception("Dimensions of b or alpha do not match the system");
  }

  /**
   * Note: T_inv and Q have only zeroes on indices reaching out the quadratic
   * size of offline.dimA(), which is the size of the original non refined grid
   */

  // stores interim values
  sgpp::base::DataVector QTQb(Q.getNrows(), 0.0);
  sgpp::base::DataVector interim(Q.getNcols(), 0.0);

  // start with Q^t * x, doing this manually until sgpp supports transpose operations
  for (size_t i = 0; i < Q.getNcols(); i++) {
    double acc = 0.0;
    for (size_t j = 0; j < Q.getNrows(); j++) {
      acc += Q.get(j, i) * b.get(j);  // note: i and j are switched in Q
    }
    QTQb.set(i, acc);
  }

  // T_inv * Q^t * x
  T_inv.mult(QTQb, interim);

  // Q * T_inv * Q^t * x
  if (B.getNcols() != 1) {
    Q.mult(interim, QTQb);
  } else if (Q.getNcols() == alpha.getSize()) {
    Q.mult(interim, alpha);
  } else {
    throw sgpp::base::data_exception("dimension of alpha doesn't match dimension of Q*T_inv*Q^t");
  }

  // if B is of size 1, then no refinement has taken place and
  // alpha is just Q * T_inv * Q^t * b
  if (B.getNcols() != 1) {
    // B * b
    B.mult(b, alpha);

    // Q*T_inv*Q^t*b + B*b
    QTQb.resize(alpha.getSize(), 0.0);
    alpha.add(QTQb);
  }
}
}  // namespace datadriven
}  // namespace sgpp
