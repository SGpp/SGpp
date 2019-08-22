// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/solver/sle/BiCGStabSP.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace solver {

BiCGStabSP::BiCGStabSP(size_t imax, float epsilon) : sgpp::solver::SLESolverSP(imax, epsilon) {}

BiCGStabSP::~BiCGStabSP() {}

void BiCGStabSP::solve(sgpp::base::OperationMatrixSP& SystemMatrix, sgpp::base::DataVectorSP& alpha,
                       sgpp::base::DataVectorSP& b, bool reuse, bool verbose, float max_threshold) {
  this->nIterations = 1;
  float epsilonSqd = this->myEpsilon * this->myEpsilon;

  if (reuse == false) {
    // Choose x0
    alpha.setAll(0.0f);
  }

  // Calculate r0
  sgpp::base::DataVectorSP r(alpha.getSize());
  SystemMatrix.mult(alpha, r);
  r.sub(b);

  float delta_0 = r.dotProduct(r) * epsilonSqd;
  float delta = 0.0f;

  if (verbose == true) {
    std::cout << "delta_0 " << delta_0 << std::endl;
  }

  // Choose r0 as r
  sgpp::base::DataVectorSP rZero(r);
  // Set p as r0
  sgpp::base::DataVectorSP p(rZero);

  float rho = rZero.dotProduct(r);
  float rho_new = 0.0f;
  float sigma = 0.0f;
  float a = 0.0f;
  float omega = 0.0f;
  float beta = 0.0f;

  sgpp::base::DataVectorSP s(alpha.getSize());
  sgpp::base::DataVectorSP v(alpha.getSize());
  sgpp::base::DataVectorSP w(alpha.getSize());

  s.setAll(0.0);
  v.setAll(0.0);
  w.setAll(0.0);

  while (this->nIterations < this->nMaxIterations) {
    // s  = Ap
    s.setAll(0.0f);
    SystemMatrix.mult(p, s);

    // std::cout << "s " << s.get(0) << " " << s.get(1)  << std::endl;

    sigma = s.dotProduct(rZero);

    if (std::abs(sigma) == 0.0f) {
      break;
    }

    a = rho / sigma;

    // w = r - a*s
    w = r;
    w.axpy((-1.0f) * a, s);

    // v = Aw
    v.setAll(0.0f);
    SystemMatrix.mult(w, v);

    // std::cout << "v " << v.get(0) << " " << v.get(1)  << std::endl;

    omega = (v.dotProduct(w)) / (v.dotProduct(v));

    // x = x - a*p - omega*w
    alpha.axpy((-1.0f) * a, p);
    alpha.axpy((-1.0f) * omega, w);

    // r = r - a*s - omega*v
    r.axpy((-1.0f) * a, s);
    r.axpy((-1.0f) * omega, v);

    rho_new = r.dotProduct(rZero);

    delta = r.dotProduct(r);

    if (verbose == true) {
      std::cout << "delta: " << delta << std::endl;
    }

    this->residuum = delta;

    // Stop in case of better accuracy
    if (delta < delta_0 || delta < max_threshold) {
      break;
    }

    beta = (rho_new / rho) * (a / omega);
    rho = rho_new;

    // p = r + beta*(p - omega*s)
    p.axpy((-1.0f) * omega, s);
    p.mult(beta);
    p.add(r);

    this->nIterations++;
  }
}

}  // namespace solver
}  // namespace sgpp
