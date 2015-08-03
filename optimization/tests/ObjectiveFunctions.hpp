#ifndef OBJECTIVE_FUNCTIONS_HPP
#define OBJECTIVE_FUNCTIONS_HPP

#include <cmath>

#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/optimization/function/ObjectiveGradient.hpp>
#include <sgpp/optimization/function/ObjectiveHessian.hpp>
#include <sgpp/optimization/function/ConstraintFunction.hpp>
#include <sgpp/optimization/function/ConstraintGradient.hpp>

class ExampleFunction : public SGPP::optimization::ObjectiveFunction {
  public:
    ExampleFunction();
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const;
};

class ExampleGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    ExampleGradient();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};

class ExampleHessian : public SGPP::optimization::ObjectiveHessian {
  public:
    ExampleHessian();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian);
    virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const;
};



class SphereGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    SphereGradient(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};

class SphereHessian : public SGPP::optimization::ObjectiveHessian {
  public:
    SphereHessian(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian);
    virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const;
};



class G3ObjectiveFunction : public SGPP::optimization::ObjectiveFunction {
  public:
    G3ObjectiveFunction(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const;
};

class G3ObjectiveGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    G3ObjectiveGradient(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};



class G3ConstraintFunction : public SGPP::optimization::ConstraintFunction {
  public:
    G3ConstraintFunction(size_t d);
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value);
};

class G3ConstraintGradient : public SGPP::optimization::ConstraintGradient {
  public:
    G3ConstraintGradient(size_t d);
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value,
              SGPP::base::DataMatrix& gradient);
};



class G8ObjectiveFunction : public SGPP::optimization::ObjectiveFunction {
  public:
    G8ObjectiveFunction();
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const;
};

class G8ObjectiveGradient : public SGPP::optimization::ObjectiveGradient {
  public:
    G8ObjectiveGradient();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};



class G8ConstraintFunction : public SGPP::optimization::ConstraintFunction {
  public:
    G8ConstraintFunction();
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value);
};

class G8ConstraintGradient : public SGPP::optimization::ConstraintGradient {
  public:
    G8ConstraintGradient();
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value,
              SGPP::base::DataMatrix& gradient);
};

#endif /* OBJECTIVE_FUNCTIONS_HPP */
