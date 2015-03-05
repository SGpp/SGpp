// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef X86_MIC_SYMMETRIC
#include <mpi.h>
#endif
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <cstdio>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    ConjugateGradients::ConjugateGradients(size_t imax, float_t epsilon) : SLESolver(imax, epsilon) {
    }

    ConjugateGradients::~ConjugateGradients() {
    }

    void ConjugateGradients::solve(SGPP::base::OperationMatrix& SystemMatrix, SGPP::base::DataVector& alpha, SGPP::base::DataVector& b, bool reuse, bool verbose, float_t max_threshold) {
      this->starting();

      if (verbose == true) {
        std::cout << "Starting Conjugated Gradients" << std::endl;
      }

      // needed for residuum calculation
      float_t epsilonSquared = this->myEpsilon * this->myEpsilon;
      // number off current iterations
      this->nIterations = 0;

      // define temporal vectors
      SGPP::base::DataVector temp(alpha.getSize());
      SGPP::base::DataVector q(alpha.getSize());
      SGPP::base::DataVector r(b);

      float_t delta_0 = 0.0;
      float_t delta_old = 0.0;
      float_t delta_new = 0.0;
      float_t beta = 0.0;
      float_t a = 0.0;

      if (verbose == true) {
        std::cout << "All temp variables used in CG have been initialized" << std::endl;
      }

      if (reuse == true) {
        q.setAll(0.0);
        SystemMatrix.mult(q, temp);
        r.sub(temp);
        delta_0 = r.dotProduct(r) * epsilonSquared;
      } else {
        alpha.setAll(0.0);
      }


      // calculate the starting residuum
      SystemMatrix.mult(alpha, temp);

      r.sub(temp);

      SGPP::base::DataVector d(r);

      delta_old = 0.0;
      delta_new = r.dotProduct(r);

      if (reuse == false) {
        delta_0 = delta_new * epsilonSquared;
      }

      this->residuum = (delta_0 / epsilonSquared);
      this->calcStarting();

      if (verbose == true) {
        std::cout << "Starting norm of residuum: " << (delta_0 / epsilonSquared) << std::endl;
        std::cout << "Target norm:               " << (delta_0) << std::endl;
      }

      while ((this->nIterations < this->nMaxIterations) && (delta_new > delta_0) && (delta_new > max_threshold)) {
        // q = A*d
        SystemMatrix.mult(d, q);

        // a = d_new / d.q
        a = delta_new / d.dotProduct(q);

        // x = x + a*d
        alpha.axpy(a, d);

        // Why ????
        if ((this->nIterations % 50) == 0) {
          // r = b - A*x
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

#ifdef X86_MIC_SYMMETRIC
        MPI_Bcast(&delta_new, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        this->residuum = delta_new;
        this->iterationComplete();

        if (verbose == true) {
          std::cout << "delta: " << delta_new << std::endl;
        }

        d.mult(beta);
        d.add(r);

        this->nIterations++;
      }

      this->residuum = delta_new;
      this->complete();

      if (verbose == true) {
        std::cout << "Number of iterations: " << this->nIterations << " (max. " << this->nMaxIterations << ")" << std::endl;
        std::cout << "Final norm of residuum: " << delta_new << std::endl;
      }
    }

    void ConjugateGradients::starting() {
    }

    void ConjugateGradients::calcStarting() {
    }

    void ConjugateGradients::iterationComplete() {
    }

    void ConjugateGradients::complete() {
    }

  }
}