// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FISTA_HPP
#define FISTA_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/solver/sle/fista/FistaBase.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>

#include <cmath>
#include <limits>
#include <type_traits>

namespace sgpp {
namespace solver {

/**
* @brief Fast Iterative Shrinkage Tresholding Algorithm is a solver for least-squares problems.
* @details
* It can solve all problems that are seperable in a least-squares part and a convex,
* not necessarily smooth other function.
* The other function is a template argument.
* FISTA is an optimal first order method for this problem class.
*/
template <typename F>
class Fista : public FistaBase {
  static_assert(std::is_base_of<RegularizationFunction, F>::value,
                "Template argument for fista is not a subtype of regularizationFunction!");

 public:
  /**
   * @brief Fista
   * @param g is the regularization function.
   */
  explicit Fista(F g) : g(g) {}
  virtual ~Fista() {}
  /**
   * @brief solve solves the problem.
   * @param op
   * @param weights is the first guess for the solution
   * @param classes is the target vector
   * @param maxIt is the maximum number of iterations
   * @param threshold is the desired accuracy
   * @param L is a guess for the Lipschitz number of the gradient of the least
   * squares part. It is used to improve the speed of the linesearch, and should
   * be 0.5 for the first iteration. Has to be positive.
   */
  void solve(base::OperationMultipleEval& op, base::DataVector& weights,
             const base::DataVector& classes, size_t maxIt, double threshold,
             double L = 0.5) override {
    // Parameters for linesearch
    // The choices made here correspond to values that work well in practice, other
    // values work as well, as long as the constraints given below still hold.
    // Larger values lead to a faster linesearch, but to worse stepsize estimates.
    const double eta = 2;  // eta > 1.0

    // Parameters for convergence check
    auto curMSE = std::numeric_limits<double>::max();
    auto priorMSE = 0.0;
    size_t curIt = 0;

    // Initial values
    base::DataVector y = weights;
    double momentum = 1;  // corresponds to t in paper

    // Variables that are reused in each iteration
    double momentumBefore = momentum;
    base::DataVector weightsBefore = base::DataVector(weights.getSize());
    auto errors = base::DataVector(classes.getSize());
    auto gradient = base::DataVector(weights.getSize());

    while (isNotConverged(curIt, maxIt, priorMSE, curMSE, threshold)) {
      priorMSE = curMSE;
      momentumBefore = momentum;
      weightsBefore = weights;

      L = L / eta;
      // Try to find the smallest L(ipschitz constant) for which our approximation works
      // First iteration is L, then L*eta.
      // Mathematically speaking, we try to find the smallest integer i \elem {1,2,...}
      // for which L*eta^i gives rise to a reasonable Lipschitz constant.
      evalErrors(op, y, classes, errors);
      evalGradient(op, errors, gradient);
      curMSE = evalResidual(errors);
      do {
        L *= eta;
        weights = evalProxGrad(y, gradient, L);              // do the step!
      } while (isNotLipschitz(weights, y, op, classes, L));  // F(prox) < Q_L(prox, weights)

      momentum = 0.5 * (1 + std::sqrt(1 + 4 * momentum * momentum));
      auto weightsDifference = weights;
      weightsDifference.sub(weightsBefore);
      const double multiplicator = (momentumBefore - 1) / momentum;
      weightsDifference.mult(multiplicator);
      y = weights;
      y.add(weightsDifference);

      ++curIt;
    }
  }

 private:
  F g;

  // Ax - b
  void evalErrors(base::OperationMultipleEval& op, base::DataVector weights,
                  const base::DataVector& b, base::DataVector& resultErrors) {
    op.mult(weights, resultErrors);
    resultErrors.sub(b);
  }

  // A^T(Ax-b)
  void evalGradient(base::OperationMultipleEval& op, base::DataVector& errors,
                    base::DataVector& resultGradient) {
    op.multTranspose(errors, resultGradient);
  }

  // 0.5 * || Ax -b ||^2
  double evalResidual(base::DataVector errors) {
    errors.sqr();
    return 0.5 * errors.sum();
  }

  // prox(weights - L^-1 gradf(a), lambda L^-1
  base::DataVector evalProxGrad(base::DataVector weights, base::DataVector gradient, double L) {
    const double stepsize = 1.0 / L;
    gradient.mult(stepsize);
    weights.sub(gradient);
    return g.prox(weights, stepsize);
  }

  // F(x) = f(x) + g(x)
  double evalGoal(double residual, const base::DataVector& weights) {
    const double gEval = g.eval(weights);
    return residual + gEval;
  }

  double evalUpperBound(base::OperationMultipleEval& op, const base::DataVector& x,
                        const base::DataVector& y, const base::DataVector& b, double L) {
    base::DataVector errors = base::DataVector(b.getSize());
    evalErrors(op, y, b, errors);
    const double fEval = evalResidual(errors);
    const double gEval = g.eval(x);
    auto gradient = base::DataVector(x.getSize());
    evalGradient(op, errors, gradient);
    auto x_min_y = x;
    x_min_y.sub(y);
    const double innerProd = gradient.dotProduct(x_min_y);
    const double x_min_y_l2 = (L / 2.0) * x_min_y.dotProduct(x_min_y);
    return fEval + innerProd + x_min_y_l2 + gEval;
  }

  bool isNotConverged(size_t curIteration, size_t maxIterations, double priorMSE, double curMSE,
                      double threshold) {
    const bool isGoodEnough = std::abs(priorMSE - curMSE) < threshold;
    const bool isFinalStep = curIteration >= maxIterations;
    return !(isGoodEnough || isFinalStep);
  }

  bool isNotLipschitz(const base::DataVector& weights, const base::DataVector& y,
                      base::OperationMultipleEval& op, const base::DataVector& b, double L) {
    base::DataVector errors = base::DataVector(b.getSize());
    evalErrors(op, weights, b, errors);
    const double residual = evalResidual(errors);
    const double goalResult = evalGoal(residual, weights);
    const double upperBoundResult = evalUpperBound(op, weights, y, b, L);
    return goalResult > upperBoundResult;
  }
};

}  // namespace solver
}  // namespace sgpp

#endif  // FISTA_HPP
