#ifndef OBJECTIVE_FUNCTIONS_HPP
#define OBJECTIVE_FUNCTIONS_HPP

#include <cmath>

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionHessian.hpp>

using namespace SGPP;
using namespace SGPP::optimization;

class ExampleFunction : public ScalarFunction {
  public:
    ExampleFunction();
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ScalarFunction>& clone) const;
};

class ExampleGradient : public ScalarFunctionGradient {
  public:
    ExampleGradient();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const;
};

class ExampleHessian : public ScalarFunctionHessian {
  public:
    ExampleHessian();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian);
    virtual void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const;
};



class SphereGradient : public ScalarFunctionGradient {
  public:
    SphereGradient(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const;
};

class SphereHessian : public ScalarFunctionHessian {
  public:
    SphereHessian(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian);
    virtual void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const;
};



class G3ObjectiveFunction : public ScalarFunction {
  public:
    G3ObjectiveFunction(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ScalarFunction>& clone) const;
};

class G3ObjectiveGradient : public ScalarFunctionGradient {
  public:
    G3ObjectiveGradient(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const;
};



class G3ConstraintFunction : public VectorFunction {
  public:
    G3ConstraintFunction(size_t d);
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value);
    virtual void clone(std::unique_ptr<VectorFunction>& clone) const;
};

class G3ConstraintGradient : public VectorFunctionGradient {
  public:
    G3ConstraintGradient(size_t d);
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value,
              SGPP::base::DataMatrix& gradient);
    virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const;
};



class G8ObjectiveFunction : public ScalarFunction {
  public:
    G8ObjectiveFunction();
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ScalarFunction>& clone) const;
};

class G8ObjectiveGradient : public ScalarFunctionGradient {
  public:
    G8ObjectiveGradient();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const;
};



class G8ConstraintFunction : public VectorFunction {
  public:
    G8ConstraintFunction();
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value);
    virtual void clone(std::unique_ptr<VectorFunction>& clone) const;
};

class G8ConstraintGradient : public VectorFunctionGradient {
  public:
    G8ConstraintGradient();
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value,
              SGPP::base::DataMatrix& gradient);
    virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const;
};

#endif /* OBJECTIVE_FUNCTIONS_HPP */
