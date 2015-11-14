#include "ObjectiveFunctions.hpp"

ExampleFunction::ExampleFunction() : ScalarFunction(2) {
}

SGPP::float_t ExampleFunction::eval(const SGPP::base::DataVector& x) {
  // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
  } else {
    return INFINITY;
  }
}

void ExampleFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(new ExampleFunction(*this));
}

ExampleGradient::ExampleGradient() : ScalarFunctionGradient(2) {
}

SGPP::float_t ExampleGradient::eval(const SGPP::base::DataVector& x,
                                    SGPP::base::DataVector& gradient) {
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    gradient[0] = 8.0 * std::cos(8.0 * x[0]);
    gradient[1] = 7.0 * std::cos(7.0 * x[1]);
    return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
  } else {
    return INFINITY;
  }
}

void ExampleGradient::clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
  clone = std::unique_ptr<ScalarFunctionGradient>(new ExampleGradient(*this));
}

ExampleHessian::ExampleHessian() : ScalarFunctionHessian(2) {
}

SGPP::float_t ExampleHessian::eval(const SGPP::base::DataVector& x,
                                   SGPP::base::DataVector& gradient,
                                   SGPP::base::DataMatrix& hessian) {
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    gradient[0] = 8.0 * std::cos(8.0 * x[0]);
    gradient[1] = 7.0 * std::cos(7.0 * x[1]);
    hessian(0, 0) = -64.0 * std::sin(8.0 * x[0]);
    hessian(0, 1) = 0.0;
    hessian(1, 0) = 0.0;
    hessian(1, 1) = -49.0 * std::sin(7.0 * x[1]);
    return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
  } else {
    return INFINITY;
  }
}

void ExampleHessian::clone(std::unique_ptr<ScalarFunctionHessian>& clone) const {
  clone = std::unique_ptr<ScalarFunctionHessian>(new ExampleHessian(*this));
}



SphereGradient::SphereGradient(size_t d) : ScalarFunctionGradient(d) {
}

SGPP::float_t SphereGradient::eval(const SGPP::base::DataVector& x,
                                   SGPP::base::DataVector& gradient) {
  SGPP::float_t result;

  for (size_t t = 0; t < d; t++) {
    if ((x[t] < 0.0) || (x[t] > 1.0)) {
      return INFINITY;
    }

    const SGPP::float_t xt = 10.0 * x[t] - 1.0;
    result += xt * xt;
    gradient[t] = 20.0 * xt;
  }

  return result;
}

void SphereGradient::clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
  clone = std::unique_ptr<ScalarFunctionGradient>(new SphereGradient(*this));
}

SphereHessian::SphereHessian(size_t d) : ScalarFunctionHessian(d) {
}

SGPP::float_t SphereHessian::eval(const SGPP::base::DataVector& x,
                                  SGPP::base::DataVector& gradient,
                                  SGPP::base::DataMatrix& hessian) {
  SGPP::float_t result;

  for (size_t t = 0; t < d; t++) {
    if ((x[t] < 0.0) || (x[t] > 1.0)) {
      return INFINITY;
    }

    const SGPP::float_t xt = 10.0 * x[t] - 1.0;

    result += xt * xt;
    gradient[t] = 20.0 * xt;
    hessian(t, t) = 20.0;

    for (size_t t2 = 0; t2 < d; t2++) {
      if (t != t2) {
        hessian(t, t2) = 0.0;
      }
    }
  }

  return result;
}

void SphereHessian::clone(std::unique_ptr<ScalarFunctionHessian>& clone) const {
  clone = std::unique_ptr<ScalarFunctionHessian>(new SphereHessian(*this));
}



DeformedLinearPhiFunction::DeformedLinearPhiFunction(size_t d) :
  VectorFunction(d, d),
  eigenvalues(d) {
  for (size_t t = 0; t < d; t++) {
    eigenvalues[t] = std::pow(10.0, t);
  }
}

void DeformedLinearPhiFunction::eval(const SGPP::base::DataVector& x,
                                     SGPP::base::DataVector& value) {
  for (size_t t = 0; t < d; t++) {
    if ((x[t] < 0.0) || (x[t] > 1.0)) {
      value.setAll(INFINITY);
      return;
    }

    value[t] = eigenvalues[t] * (x[t] - 0.1);
  }
}

void DeformedLinearPhiFunction::clone(
  std::unique_ptr<VectorFunction>& clone) const {
  clone = std::unique_ptr<VectorFunction>(
            new DeformedLinearPhiFunction(*this));
}



DeformedLinearPhiGradient::DeformedLinearPhiGradient(size_t d) :
  VectorFunctionGradient(d, d),
  eigenvalues(d) {
  for (size_t t = 0; t < d; t++) {
    eigenvalues[t] = std::pow(10.0, t);
  }
}

void DeformedLinearPhiGradient::eval(const SGPP::base::DataVector& x,
                                     SGPP::base::DataVector& value,
                                     SGPP::base::DataMatrix& gradient) {
  for (size_t t = 0; t < d; t++) {
    if ((x[t] < 0.0) || (x[t] > 1.0)) {
      value.setAll(INFINITY);
      gradient.setAll(NAN);
      return;
    }

    value[t] = eigenvalues[t] * (x[t] - 0.1);

    for (size_t t2 = 0; t2 < d; t2++) {
      gradient(t, t2) = ((t == t2) ? eigenvalues[t] : 0.0);
    }
  }
}

void DeformedLinearPhiGradient::clone(
  std::unique_ptr<VectorFunctionGradient>& clone) const {
  clone = std::unique_ptr<VectorFunctionGradient>(
            new DeformedLinearPhiGradient(*this));
}



G3ObjectiveFunction::G3ObjectiveFunction(size_t d) : ScalarFunction(d) {
}

SGPP::float_t G3ObjectiveFunction::eval(const SGPP::base::DataVector& x) {
  const SGPP::float_t dDbl = static_cast<SGPP::float_t>(d);
  SGPP::float_t fx = std::pow(dDbl, dDbl / 2.0);

  for (size_t t = 0; t < d; t++) {
    if ((x[0] >= 0.0) && (x[1] <= 1.0)) {
      fx *= x[t];
    } else {
      return INFINITY;
    }
  }

  return fx;
}

void G3ObjectiveFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(
            new G3ObjectiveFunction(*this));
}

G3ObjectiveGradient::G3ObjectiveGradient(size_t d) : ScalarFunctionGradient(d) {
}

SGPP::float_t G3ObjectiveGradient::eval(const SGPP::base::DataVector& x,
                                        SGPP::base::DataVector& gradient) {
  const SGPP::float_t dDbl = static_cast<SGPP::float_t>(d);
  SGPP::float_t fx = std::pow(dDbl, dDbl / 2.0);

  for (size_t t = 0; t < d; t++) {
    if ((x[0] >= 0.0) && (x[1] <= 1.0)) {
      gradient[t] = fx;
    } else {
      return INFINITY;
    }
  }

  for (size_t t = 0; t < d; t++) {
    for (size_t t2 = 0; t2 < d; t2++) {
      if (t2 != t) {
        gradient[t2] *= x[t];
      } else {
        fx *= x[t];
      }
    }
  }

  return fx;
}

void G3ObjectiveGradient::clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
  clone = std::unique_ptr<ScalarFunctionGradient>(
            new G3ObjectiveGradient(*this));
}



G3ConstraintFunction::G3ConstraintFunction(size_t d) : VectorFunction(d, 1) {
}

void G3ConstraintFunction::eval(const SGPP::base::DataVector& x,
                                SGPP::base::DataVector& value) {
  SGPP::float_t gx = -1.0;

  for (size_t t = 0; t < d; t++) {
    if ((x[0] >= 0.0) && (x[1] <= 1.0)) {
      gx += x[t] * x[t];
    } else {
      value[0] = INFINITY;
      return;
    }
  }

  value[0] = gx;
}

void G3ConstraintFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::unique_ptr<VectorFunction>(
            new G3ConstraintFunction(*this));
}

G3ConstraintGradient::G3ConstraintGradient(size_t d) : VectorFunctionGradient(d, 1) {
}

void G3ConstraintGradient::eval(const SGPP::base::DataVector& x,
                                SGPP::base::DataVector& value,
                                SGPP::base::DataMatrix& gradient) {
  SGPP::float_t gx = -1.0;

  for (size_t t = 0; t < d; t++) {
    if ((x[0] >= 0.0) && (x[1] <= 1.0)) {
      gx += x[t] * x[t];
      gradient(0, t) = 2.0 * x[t];
    } else {
      value[0] = INFINITY;
      return;
    }
  }

  value[0] = gx;
}

void G3ConstraintGradient::clone(std::unique_ptr<VectorFunctionGradient>& clone) const {
  clone = std::unique_ptr<VectorFunctionGradient>(
            new G3ConstraintGradient(*this));
}



G8ObjectiveFunction::G8ObjectiveFunction() : ScalarFunction(2) {
}

SGPP::float_t G8ObjectiveFunction::eval(const SGPP::base::DataVector& x) {
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    const SGPP::float_t x0 = 10.0 * x[0];
    const SGPP::float_t x1 = 10.0 * x[1];
    const SGPP::float_t fx = -std::pow(std::sin(2.0 * M_PI * x0), 3.0) *
                             std::sin(2.0 * M_PI * x1) / (std::pow(x0, 3.0) * (x0 + x1));
    return fx;
  } else {
    return INFINITY;
  }
}

void G8ObjectiveFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(
            new G8ObjectiveFunction(*this));
}

G8ObjectiveGradient::G8ObjectiveGradient() : ScalarFunctionGradient(2) {
}

SGPP::float_t G8ObjectiveGradient::eval(const SGPP::base::DataVector& x,
                                        SGPP::base::DataVector& gradient) {
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    const SGPP::float_t x0 = 10.0 * x[0];
    const SGPP::float_t x1 = 10.0 * x[1];
    const SGPP::float_t fx = -std::pow(std::sin(2.0 * M_PI * x0), 3.0) *
                             std::sin(2.0 * M_PI * x1) / (std::pow(x0, 3.0) * (x0 + x1));
    gradient[0] = 6.0 * M_PI * std::cos(2.0 * M_PI * x0) *
                  std::pow(std::sin(2.0 * M_PI * x0), 2.0) *
                  std::sin(2.0 * M_PI * x1) / ((x0 + x1) * std::pow(x0, 3.0)) -
                  3.0 * std::pow(std::sin(2.0 * M_PI * x0), 3.0) *
                  std::sin(2.0 * M_PI * x1) / ((x0 + x1) * std::pow(x0, 4.0)) -
                  std::pow(std::sin(2.0 * M_PI * x0), 3.0) *
                  std::sin(2.0 * M_PI * x1) /
                  (std::pow(x0 + x1, 2.0) * std::pow(x0, 3.0));
    gradient[1] = 2.0 * M_PI * std::cos(2.0 * M_PI * x1) *
                  std::pow(std::sin(2.0 * M_PI * x0), 3.0) /
                  ((x0 + x1) * std::pow(x0, 3.0)) -
                  std::pow(std::sin(2.0 * M_PI * x0), 3.0) *
                  std::sin(2.0 * M_PI * x1) /
                  (std::pow(x0 + x1, 2.0) * std::pow(x0, 3.0));
    gradient[0] *= -10.0;
    gradient[1] *= -10.0;
    return fx;
  } else {
    return INFINITY;
  }
}

void G8ObjectiveGradient::clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
  clone = std::unique_ptr<ScalarFunctionGradient>(
            new G8ObjectiveGradient(*this));
}



G8ConstraintFunction::G8ConstraintFunction() : VectorFunction(2, 2) {
}

void G8ConstraintFunction::eval(const SGPP::base::DataVector& x,
                                SGPP::base::DataVector& value) {
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    const SGPP::float_t x0 = 10.0 * x[0];
    const SGPP::float_t x1 = 10.0 * x[1];
    value[0] = std::pow(x0, 2.0) - x1 + 1.0;
    value[1] = 1.0 - x0 + std::pow(x1 - 4.0, 2.0);
  } else {
    value[0] = INFINITY;
    value[1] = INFINITY;
    return;
  }
}

void G8ConstraintFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::unique_ptr<VectorFunction>(
            new G8ConstraintFunction(*this));
}

G8ConstraintGradient::G8ConstraintGradient() : VectorFunctionGradient(2, 2) {
}

void G8ConstraintGradient::eval(const SGPP::base::DataVector& x,
                                SGPP::base::DataVector& value,
                                SGPP::base::DataMatrix& gradient) {
  if ((x[0] >= 0.0) && (x[0] <= 1.0) &&
      (x[1] >= 0.0) && (x[1] <= 1.0)) {
    const SGPP::float_t x0 = 10.0 * x[0];
    const SGPP::float_t x1 = 10.0 * x[1];
    value[0] = std::pow(x0, 2.0) - x1 + 1.0;
    value[1] = 1.0 - x0 + std::pow(x1 - 4.0, 2.0);
    gradient(0, 0) = 2.0 * x0 * 10.0;
    gradient(0, 1) = -1.0 * 10.0;
    gradient(1, 0) = -1.0 * 10.0;
    gradient(1, 1) = 2.0 * (x1 - 4.0) * 10.0;
  } else {
    value[0] = INFINITY;
    value[1] = INFINITY;
    return;
  }
}

void G8ConstraintGradient::clone(std::unique_ptr<VectorFunctionGradient>& clone) const {
  clone = std::unique_ptr<VectorFunctionGradient>(
            new G8ConstraintGradient(*this));
}
