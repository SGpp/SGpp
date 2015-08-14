#ifndef OBJECTIVE_FUNCTIONS_HPP
#define OBJECTIVE_FUNCTIONS_HPP

#include <cmath>

#include <sgpp/optimization/function/scalar/ObjectiveFunction.hpp>
#include <sgpp/optimization/function/vector/ConstraintFunction.hpp>

using namespace SGPP;
using namespace SGPP::optimization;

class ExampleFunction : public ObjectiveFunction {
  public:
    ExampleFunction();
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const;
};

class ExampleGradient : public ObjectiveGradient {
  public:
    ExampleGradient();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};

class ExampleHessian : public ObjectiveHessian {
  public:
    ExampleHessian();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian);
    virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const;
};



class SphereGradient : public ObjectiveGradient {
  public:
    SphereGradient(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};

class SphereHessian : public ObjectiveHessian {
  public:
    SphereHessian(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient,
                       SGPP::base::DataMatrix& hessian);
    virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const;
};



class G3ObjectiveFunction : public ObjectiveFunction {
  public:
    G3ObjectiveFunction(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const;
};

class G3ObjectiveGradient : public ObjectiveGradient {
  public:
    G3ObjectiveGradient(size_t d);
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};



class G3ConstraintFunction : public ConstraintFunction {
  public:
    G3ConstraintFunction(size_t d);
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value);
    virtual void clone(std::unique_ptr<ConstraintFunction>& clone) const;
};

class G3ConstraintGradient : public ConstraintGradient {
  public:
    G3ConstraintGradient(size_t d);
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value,
              SGPP::base::DataMatrix& gradient);
    virtual void clone(std::unique_ptr<ConstraintGradient>& clone) const;
};



class G8ObjectiveFunction : public ObjectiveFunction {
  public:
    G8ObjectiveFunction();
    SGPP::float_t eval(const SGPP::base::DataVector& x);
    virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const;
};

class G8ObjectiveGradient : public ObjectiveGradient {
  public:
    G8ObjectiveGradient();
    SGPP::float_t eval(const SGPP::base::DataVector& x,
                       SGPP::base::DataVector& gradient);
    virtual void clone(std::unique_ptr<ObjectiveGradient>& clone) const;
};



class G8ConstraintFunction : public ConstraintFunction {
  public:
    G8ConstraintFunction();
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value);
    virtual void clone(std::unique_ptr<ConstraintFunction>& clone) const;
};

class G8ConstraintGradient : public ConstraintGradient {
  public:
    G8ConstraintGradient();
    void eval(const SGPP::base::DataVector& x,
              SGPP::base::DataVector& value,
              SGPP::base::DataMatrix& gradient);
    virtual void clone(std::unique_ptr<ConstraintGradient>& clone) const;
};

#endif /* OBJECTIVE_FUNCTIONS_HPP */
