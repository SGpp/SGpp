// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>

namespace SGPP {
namespace solver {

BiCGStab::BiCGStab(size_t imax, float_t epsilon) : SLESolver(imax, epsilon) {}

BiCGStab::~BiCGStab() {}

void BiCGStab::solve(SGPP::base::OperationMatrix& SystemMatrix, SGPP::base::DataVector& alpha,
                     SGPP::base::DataVector& b, bool reuse, bool verbose, float_t max_threshold) {
  this->nIterations = 1;
  float_t epsilonSqd = this->myEpsilon * this->myEpsilon;

  if (reuse == false) {
    // Choose x0
    alpha.setAll(0.0);
  }

  // Calculate r0
  SGPP::base::DataVector r(alpha.getSize());
  SystemMatrix.mult(alpha, r);
  r.sub(b);

  float_t delta_0 = r.dotProduct(r) * epsilonSqd;
  float_t delta = 0.0;

  if (verbose == true) {
    std::cout << "delta_0 " << delta_0 << std::endl;
  }

  // Choose r0 as r
  SGPP::base::DataVector rZero(r);
  // Set p as r0
  SGPP::base::DataVector p(rZero);

  float_t rho = rZero.dotProduct(r);
  float_t rho_new = 0.0;
  float_t sigma = 0.0;
  float_t a = 0.0;
  float_t omega = 0.0;
  float_t beta = 0.0;

  SGPP::base::DataVector s(alpha.getSize());
  SGPP::base::DataVector v(alpha.getSize());
  SGPP::base::DataVector w(alpha.getSize());

  s.setAll(0.0);
  v.setAll(0.0);
  w.setAll(0.0);

  while (this->nIterations < this->nMaxIterations) {
    // s  = Ap
    s.setAll(0.0);
    SystemMatrix.mult(p, s);

    // std::cout << "s " << s.get(0) << " " << s.get(1)  << std::endl;

    sigma = s.dotProduct(rZero);

    if (fabs(sigma) == 0.0) {
      break;
    }

    a = rho / sigma;

    // w = r - a*s
    w = r;
    w.axpy((-1.0) * a, s);

    // v = Aw
    v.setAll(0.0);
    SystemMatrix.mult(w, v);

    // std::cout << "v " << v.get(0) << " " << v.get(1)  << std::endl;

    omega = (v.dotProduct(w)) / (v.dotProduct(v));

    // x = x - a*p - omega*w
    alpha.axpy((-1.0) * a, p);
    alpha.axpy((-1.0) * omega, w);

    // r = r - a*s - omega*v
    r.axpy((-1.0) * a, s);
    r.axpy((-1.0) * omega, v);

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
    p.axpy((-1.0) * omega, s);
    p.mult(beta);
    p.add(r);

    this->nIterations++;
  }
}

}  // namespace solver
}  // namespace SGPP
