// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#endif

#include <sgpp/parallel/solver/sle/ConjugateGradientsMPI.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace parallel {

ConjugateGradientsMPI::ConjugateGradientsMPI(size_t imax, double epsilon)
    : SGPP::solver::SLESolver(imax, epsilon) {}

ConjugateGradientsMPI::~ConjugateGradientsMPI() {}

void ConjugateGradientsMPI::solve(SGPP::base::OperationMatrix& SystemMatrix,
                                  SGPP::base::DataVector& alpha, SGPP::base::DataVector& b,
                                  bool reuse, bool verbose, double max_threshold) {
  myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);
  myGlobalMPIComm->broadcastGridCoefficientsFromRank0(b);

  //  if (myGlobalMPIComm->getMyRank() != 0)
  //  {
  //    this->waitForTask(SystemMatrix, alpha);
  //  }
  //  else
  //  {
  if (verbose == true) {
    std::cout << "Starting Conjugated Gradients" << std::endl;
  }

  //    char ctrl = 'M';
  // needed for residuum calculation
  double epsilonSquared = this->myEpsilon * this->myEpsilon;
  // number off current iterations
  this->nIterations = 0;

  // define temporal vectors
  SGPP::base::DataVector temp(alpha.getSize());
  SGPP::base::DataVector q(alpha.getSize());
  SGPP::base::DataVector r(b);

  double delta_0 = 0.0;
  double delta_old = 0.0;
  double delta_new = 0.0;
  double beta = 0.0;
  double a = 0.0;

  if (verbose == true) {
    std::cout << "All temp variables used in CG have been initialized" << std::endl;
  }

  if (reuse == true) {
    q.setAll(0.0);

    //      ctrl = 'M';
    //      myGlobalMPIComm->broadcastControlFromRank0(&ctrl);
    SystemMatrix.mult(q, temp);

    r.sub(temp);
    delta_0 = r.dotProduct(r) * epsilonSquared;
  } else {
    alpha.setAll(0.0);
  }

  // calculate the starting residuum
  //    ctrl = 'M';
  //    myGlobalMPIComm->broadcastControlFromRank0(&ctrl);
  SystemMatrix.mult(alpha, temp);

  r.sub(temp);

  SGPP::base::DataVector d(r);

  delta_old = 0.0;
  delta_new = r.dotProduct(r);

  if (reuse == false) {
    delta_0 = delta_new * epsilonSquared;
  }

  this->residuum = (delta_0 / epsilonSquared);

  if (verbose == true) {
    std::cout << "Starting norm of residuum: " << (delta_0 / epsilonSquared) << std::endl;
    std::cout << "Target norm:               " << (delta_0) << std::endl;
  }

  while ((this->nIterations < this->nMaxIterations) && (delta_new > delta_0) &&
         (delta_new > max_threshold)) {
    // q = A*d
    //      ctrl = 'M';
    //      myGlobalMPIComm->broadcastControlFromRank0(&ctrl);
    SystemMatrix.mult(d, q);

    // a = d_new / d.q
    a = delta_new / d.dotProduct(q);

    // x = x + a*d
    alpha.axpy(a, d);

    // Why ????
    if ((this->nIterations % 50) == 0) {
      // r = b - A*x
      //        ctrl = 'M';
      //        myGlobalMPIComm->broadcastControlFromRank0(&ctrl);
      SystemMatrix.mult(alpha, temp);

      r.copyFrom(b);
      r.sub(temp);
    } else {
      // r = r - a*q
      r.axpy(-a, q);
    }

    // calculate new deltas and determine beta
    delta_old = delta_new;
    delta_new = r.dotProduct(r);
    beta = delta_new / delta_old;

    this->residuum = delta_new;

    if (verbose == true) {
      std::cout << "delta: " << delta_new << std::endl;
    }

    d.mult(beta);
    d.add(r);

    this->nIterations++;
  }

  this->residuum = delta_new;

  //    ctrl = 'T';
  //    myGlobalMPIComm->broadcastControlFromRank0(&ctrl);

  if (verbose == true) {
    std::cout << "Number of iterations: " << this->nIterations << " (max. " << this->nMaxIterations
              << ")" << std::endl;
    std::cout << "Final norm of residuum: " << delta_new << std::endl;
  }

  //  }
}

void ConjugateGradientsMPI::waitForTask(SGPP::base::OperationMatrix& SystemMatrix,
                                        SGPP::base::DataVector& alpha) {
  char ctrl;
  SGPP::base::DataVector result(alpha.getSize());

  do {
    myGlobalMPIComm->broadcastControlFromRank0(&ctrl);

    if (ctrl == 'M') {
      SystemMatrix.mult(alpha, result);
    }
  } while (ctrl != 'T');
}
}  // namespace parallel
}  // namespace SGPP
