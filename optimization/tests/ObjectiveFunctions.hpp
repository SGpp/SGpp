#include <cmath>

#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/optimization/function/ObjectiveGradient.hpp>
#include <sgpp/optimization/function/ObjectiveHessian.hpp>
#include <sgpp/optimization/function/ConstraintFunction.hpp>
#include <sgpp/optimization/function/ConstraintGradient.hpp>

class ExampleFunction : public SGPP::optimization::ObjectiveFunction {
  public:
    ExampleFunction() : ObjectiveFunction(2) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x) {
      // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
      if ((x.get(0) >= 0.0) && (x.get(0) <= 1.0) &&
          (x.get(1) >= 0.0) && (x.get(1) <= 1.0)) {
        return std::sin(8.0 * x.get(0)) + std::sin(7.0 * x.get(1));
      } else {
        return INFINITY;
      }
    }

    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
      clone = std::unique_ptr<ObjectiveFunction>(new ExampleFunction(*this));
    }
};

class ExampleGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    ExampleGradient() : ObjectiveGradient(2) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient) {
      if ((x.get(0) >= 0.0) && (x.get(0) <= 1.0) &&
          (x.get(1) >= 0.0) && (x.get(1) <= 1.0)) {
        gradient[0] = 8.0 * std::cos(8.0 * x.get(0));
        gradient[1] = 7.0 * std::cos(7.0 * x.get(1));
        return std::sin(8.0 * x.get(0)) + std::sin(7.0 * x.get(1));
      } else {
        return INFINITY;
      }
    }

    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const {
      clone = std::unique_ptr<ObjectiveGradient>(new ExampleGradient(*this));
    }
};

class ExampleHessian : public SGPP::optimization::ObjectiveHessian {
  public:
    ExampleHessian() : ObjectiveHessian(2) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian) {
      if ((x.get(0) >= 0.0) && (x.get(0) <= 1.0) &&
          (x.get(1) >= 0.0) && (x.get(1) <= 1.0)) {
        gradient[0] = 8.0 * std::cos(8.0 * x.get(0));
        gradient[1] = 7.0 * std::cos(7.0 * x.get(1));
        hessian.set(0, 0, -64.0 * std::sin(8.0 * x.get(0)));
        hessian.set(0, 1, 0.0);
        hessian.set(1, 0, 0.0);
        hessian.set(1, 1, -49.0 * std::sin(7.0 * x.get(1)));
        return std::sin(8.0 * x.get(0)) + std::sin(7.0 * x.get(1));
      } else {
        return INFINITY;
      }
    }

    virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const {
      clone = std::unique_ptr<ObjectiveHessian>(new ExampleHessian(*this));
    }
};



class SphereGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    SphereGradient(size_t d) : ObjectiveGradient(d) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient) {
      SGPP::float_t result;

      for (size_t t = 0; t < d; t++) {
        if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
          return INFINITY;
        }

        const SGPP::float_t xt = 10.0 * x.get(t) - 1.0;
        result += xt * xt;
        gradient[t] = 20.0 * xt;
      }

      return result;
    }

    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const {
      clone = std::unique_ptr<ObjectiveGradient>(new SphereGradient(*this));
    }
};

class SphereHessian : public SGPP::optimization::ObjectiveHessian {
  public:
    SphereHessian(size_t d) : ObjectiveHessian(d) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian) {
      SGPP::float_t result;

      for (size_t t = 0; t < d; t++) {
        if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
          return INFINITY;
        }

        const SGPP::float_t xt = 10.0 * x.get(t) - 1.0;

        result += xt * xt;
        gradient[t] = 20.0 * xt;
        hessian.set(t, t, 20.0);

        for (size_t t2 = 0; t2 < d; t2++) {
          if (t != t2) {
            hessian.set(t, t2, 0.0);
          }
        }
      }

      return result;
    }

    virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const {
      clone = std::unique_ptr<ObjectiveHessian>(new SphereHessian(*this));
    }
};



/*class G3ObjectiveFunction : public SGPP::optimization::ObjectiveFunction {
  public:
    G3ObjectiveFunction(size_t d) : G3ObjectiveFunction(d) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x) {
      SGPP::float_t fx = std::pow(d, d / 2.0);

      for (size_t t = 0; t < d; t++) {
        if ((x.get(0) >= 0.0) && (x.get(1) <= 1.0)) {
          fx *= x.get(t);
        } else {
          return INFINITY;
        }
      }

      return fx;
    }

    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
      clone = std::unique_ptr<ObjectiveFunction>(
          new G3ObjectiveFunction(*this));
    }
};

class G3ObjectiveGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    G3ObjectiveGradient(size_t d) : ObjectiveGradient(d) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient) {
      SGPP::float_t fx = std::pow(d, d / 2.0);

      for (size_t t = 0; t < d; t++) {
        if ((x.get(0) >= 0.0) && (x.get(1) <= 1.0)) {
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

    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const {
      clone = std::unique_ptr<ObjectiveGradient>(
          new G3ObjectiveGradient(*this));
    }
};



class G3ConstraintFunction : public SGPP::optimization::ConstraintFunction {
  public:
    G3ConstraintFunction(size_t d) : ConstraintFunction(d) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x) {
      SGPP::float_t fx = std::pow(d, d / 2.0);

      for (size_t t = 0; t < d; t++) {
        if ((x.get(0) >= 0.0) && (x.get(1) <= 1.0)) {
          fx *= x.get(t);
        } else {
          return INFINITY;
        }
      }

      return fx;
    }

    virtual void clone(std::unique_ptr<ConstraintFunction>& clone) const {
      clone = std::unique_ptr<ConstraintFunction>(
          new G3ConstraintFunction(*this));
    }
};

class G3ConstraintGradient : public SGPP::optimization::ConstraintGradient {
  public:
    G3ConstraintGradient(size_t d) : ConstraintGradient(d) {
    }

    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient) {
      SGPP::float_t fx = std::pow(d, d / 2.0);

      for (size_t t = 0; t < d; t++) {
        if ((x.get(0) >= 0.0) && (x.get(1) <= 1.0)) {
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

    virtual void clone(std::unique_ptr<ConstraintGradient>& clone) const {
      clone = std::unique_ptr<ConstraintGradient>(
          new G3ConstraintGradient(*this));
    }
};*/
