// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp>

#include <iomanip>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Struct for a LagrangePolynom, used to eval it
 */
struct LagrangePolynom {
  std::vector<double> points;
  size_t point;

  double evaluate(double x) {
    double result = 1.0;
    for (size_t i = 0; i < points.size(); ++i) {
      if (i != point) {
        // TODO(holzmudd): precalculate denominator?
        result *= (x - points[i]) / (points[point] - points[i]);
      }
    }
    return result;
  }
};

/*
 *	Implementation of Gauss-Legendre quadrature
 */
class GaussLegendreQuadrature {
 public:
  explicit GaussLegendreQuadrature(size_t eDEGREE)
      : s_LegendrePolynomial(eDEGREE), eDEGREE(eDEGREE) {}

  /*! Compute the integral of a functor
   *
   *   @param a    lower limit of integration
   *   @param b    upper limit of integration
   *   @param f    the polynomial to integrate
   */
  double integrate(double a, double b, LagrangePolynom f) {
    double p = (b - a) / 2;
    double q = (b + a) / 2;
    const LegendrePolynomial& legpoly = s_LegendrePolynomial;

    double sum = 0;
    for (size_t i = 1; i <= eDEGREE; ++i) {  // TODO(holzmudd): warum hier ab 1?
      sum += legpoly.weight(i) * f.evaluate(p * legpoly.root(i) + q);
    }

    return p * sum;
  }

 private:
  /*
   * Implementation of the Legendre polynomials that form
   * the basis of this quadrature
   */
  class LegendrePolynomial {
    size_t eDEGREE;

   public:
    explicit LegendrePolynomial(size_t eDEGREE) : eDEGREE(eDEGREE) {
      _r = std::vector<double>(
          eDEGREE + 1,
          0.0);  // TODO(holzmudd): die beiden waren (eDEGREE, 0.0), wie rum ist es richtig?
      _w = std::vector<double>(eDEGREE + 1, 0.0);
      // Solve roots and weights
      for (size_t i = 0; i <= eDEGREE;
           ++i) {  // TODO(holzmudd): war i <= eDEGREE, wie rum ist es richtig?
        double dr = 1;

        // Find zero
        Evaluation eval(
            cos(M_PI * (static_cast<double>(i) - 0.25) / (static_cast<double>(eDEGREE) + 0.5)),
            eDEGREE);
        do {
          dr = eval.v() / eval.d();
          eval.evaluate(eval.x() - dr);
        } while (fabs(dr) > 2e-16);

        this->_r[i] = eval.x();
        this->_w[i] = 2 / ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
      }
    }

    double root(size_t i) const { return this->_r[i]; }
    double weight(size_t i) const { return this->_w[i]; }

   private:
    std::vector<double> _r;
    std::vector<double> _w;

    /*
     * Evaluate the value *and* derivative of the
     * Legendre polynomial
     */
    class Evaluation {
      size_t eDEGREE;

     public:
      explicit Evaluation(double x, size_t eDEGREE) : eDEGREE(eDEGREE), _x(x), _v(1), _d(0) {
        this->evaluate(x);
      }

      void evaluate(double x) {
        this->_x = x;

        double vsub1 = x;
        double vsub2 = 1;
        double f = 1 / (x * x - 1);

        for (size_t i = 2; i <= eDEGREE; ++i) {  // TODO(holzmudd): wirklich hier ab 2?
          this->_v = ((2 * static_cast<double>(i) - 1) * x * vsub1 -
                      (static_cast<double>(i) - 1) * vsub2) /
                     static_cast<double>(i);
          this->_d = static_cast<double>(i) * f * (x * this->_v - vsub1);

          vsub2 = vsub1;
          vsub1 = this->_v;
        }
      }

      double v() const { return this->_v; }
      double d() const { return this->_d; }
      double x() const { return this->_x; }

     private:
      double _x;
      double _v;
      double _d;
    };
  };

  /*
   * Pre-compute the weights and abscissae of the Legendre polynomials
   */
  LegendrePolynomial s_LegendrePolynomial;
  size_t eDEGREE;
};

double gauss(LagrangePolynom polynom) {
  size_t size = polynom.points.size() + 1;
  GaussLegendreQuadrature quad(
      size);  // TODO(holzmudd): kann man optimieren, indem man das nicht immer neu erstellt
  return quad.integrate(0.0, 1.0, polynom);
}

/**
 * Integrates the polynom
 */
double integrate(LagrangePolynom polynom) {
  // return romberg2(polynom);
  return gauss(polynom);
}

/**
 * Calculates the weight for the specific point
 */
double getWeight(std::vector<double>& points, size_t point) {
  LagrangePolynom p;
  p.points = points;
  p.point = point;

  return integrate(p);
}

/**
 * This Function calculates the weights of the given points, each weight is calculated individually
 * @param points The vector with the points, they dont need to have a specific order
 * @param weights The weights will be added to the back of this vector in the order of the points in
 * the vector with the points,
 * it is recommended to clear the weight vector before calling this function to ensure that the
 * weights are at the same position
 * as their points
 */
void calculateWeights(std::vector<double>& points, std::vector<FloatScalarVector>& weights) {
  // calc weight for each point
  for (size_t i = 0; i < points.size(); ++i) {
    weights.push_back(FloatScalarVector(getWeight(points, i)));
  }
}

QuadratureEvaluator::~QuadratureEvaluator() {}

bool QuadratureEvaluator::needsOrderedPoints() { return false; }

bool QuadratureEvaluator::needsParameter() { return false; }

void QuadratureEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  weights.clear();
  calculateWeights(xValues, weights);

  if (normalizeWeights) {
    double sum = 0.0;

    // multiply the weights with the weight function
    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(weight_function(xValues[i]));
      sum += weights[i].getValue();
    }

    double sumInv = 1.0 / sum;

    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(sumInv);
    }
  } else {
    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(weight_function(xValues[i]));
    }
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > QuadratureEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new QuadratureEvaluator(*this));
}

QuadratureEvaluator::QuadratureEvaluator(sgpp::combigrid::SingleFunction weight_function,
                                         bool normalizeWeights)
    : weight_function(weight_function), normalizeWeights(normalizeWeights) {}

QuadratureEvaluator::QuadratureEvaluator(QuadratureEvaluator const& other)
    : xValues(other.xValues),
      weights(other.weights),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights) {}

void QuadratureEvaluator::setParameter(const FloatScalarVector& param) { return; }

} /* namespace combigrid */
} /* namespace sgpp*/
