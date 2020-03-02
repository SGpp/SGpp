// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace solver {

BiCGStab::BiCGStab(size_t imax, double epsilon) : SLESolver(imax, epsilon) {}

BiCGStab::~BiCGStab() {}

void BiCGStab::solve(sgpp::base::OperationMatrix& SystemMatrix, sgpp::base::DataVector& alpha,
                     sgpp::base::DataVector& b, bool reuse, bool verbose, double max_threshold) {
  this->nIterations = 1;
  double epsilonSqd = this->myEpsilon * this->myEpsilon;

  if (reuse == false) {
    // Choose x0
    alpha.setAll(0.0);
  }

  // Calculate r0
  sgpp::base::DataVector r(alpha.getSize());
  SystemMatrix.mult(alpha, r);
  r.sub(b);

  double delta_0 = r.dotProduct(r) * epsilonSqd;
  double delta = 0.0;

  if (verbose == true) {
    std::cout << "delta_0 " << delta_0 << std::endl;
  }

  // Choose r0 as r
  sgpp::base::DataVector rZero(r);
  // Set p as r0
  sgpp::base::DataVector p(rZero);

  double rho = rZero.dotProduct(r);
  double rho_new = 0.0;
  double sigma = 0.0;
  double a = 0.0;
  double omega = 0.0;
  double beta = 0.0;

  sgpp::base::DataVector s(alpha.getSize());
  sgpp::base::DataVector v(alpha.getSize());
  sgpp::base::DataVector w(alpha.getSize());

  while (this->nIterations < this->nMaxIterations) {
    // s  = Ap
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
}  // namespace sgpp
