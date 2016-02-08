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
  virtual ~ExampleFunction() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x) override;
  virtual void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

class ExampleGradient : public ScalarFunctionGradient {
 public:
  ExampleGradient();
  virtual ~ExampleGradient() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x,
                             SGPP::base::DataVector& gradient) override;
  virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};

class ExampleHessian : public ScalarFunctionHessian {
 public:
  ExampleHessian();
  virtual ~ExampleHessian() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x,
                             SGPP::base::DataVector& gradient,
                             SGPP::base::DataMatrix& hessian) override;
  virtual void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const
  override;
};



class SphereGradient : public ScalarFunctionGradient {
 public:
  SphereGradient(size_t d);
  virtual ~SphereGradient() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x,
                             SGPP::base::DataVector& gradient) override;
  virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};

class SphereHessian : public ScalarFunctionHessian {
 public:
  SphereHessian(size_t d);
  virtual ~SphereHessian() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x,
                             SGPP::base::DataVector& gradient,
                             SGPP::base::DataMatrix& hessian) override;
  virtual void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const
  override;
};



class DeformedLinearPhiFunction : public VectorFunction {
 public:
  DeformedLinearPhiFunction(size_t d);
  virtual ~DeformedLinearPhiFunction() override;
  virtual void eval(const SGPP::base::DataVector& x,
                    SGPP::base::DataVector& value) override;
  virtual void clone(std::unique_ptr<VectorFunction>& clone) const override;

 protected:
  SGPP::base::DataVector eigenvalues;
};

class DeformedLinearPhiGradient : public VectorFunctionGradient {
 public:
  DeformedLinearPhiGradient(size_t d);
  virtual ~DeformedLinearPhiGradient() override;
  virtual void eval(const SGPP::base::DataVector& x,
                    SGPP::base::DataVector& value,
                    SGPP::base::DataMatrix& gradient) override;
  virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override;

 protected:
  SGPP::base::DataVector eigenvalues;
};



class G3ObjectiveFunction : public ScalarFunction {
 public:
  G3ObjectiveFunction(size_t d);
  virtual ~G3ObjectiveFunction() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x) override;
  virtual void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

class G3ObjectiveGradient : public ScalarFunctionGradient {
 public:
  G3ObjectiveGradient(size_t d);
  virtual ~G3ObjectiveGradient() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x,
                             SGPP::base::DataVector& gradient) override;
  virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};



class G3ConstraintFunction : public VectorFunction {
 public:
  G3ConstraintFunction(size_t d);
  virtual ~G3ConstraintFunction() override;
  virtual void eval(const SGPP::base::DataVector& x,
                    SGPP::base::DataVector& value) override;
  virtual void clone(std::unique_ptr<VectorFunction>& clone) const override;
};

class G3ConstraintGradient : public VectorFunctionGradient {
 public:
  G3ConstraintGradient(size_t d);
  virtual ~G3ConstraintGradient() override;
  virtual void eval(const SGPP::base::DataVector& x,
                    SGPP::base::DataVector& value,
                    SGPP::base::DataMatrix& gradient) override;
  virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override;
};



class G8ObjectiveFunction : public ScalarFunction {
 public:
  G8ObjectiveFunction();
  virtual ~G8ObjectiveFunction() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x) override;
  virtual void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

class G8ObjectiveGradient : public ScalarFunctionGradient {
 public:
  G8ObjectiveGradient();
  virtual ~G8ObjectiveGradient() override;
  virtual SGPP::float_t eval(const SGPP::base::DataVector& x,
                             SGPP::base::DataVector& gradient) override;
  virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};



class G8ConstraintFunction : public VectorFunction {
 public:
  G8ConstraintFunction();
  virtual ~G8ConstraintFunction() override;
  virtual void eval(const SGPP::base::DataVector& x,
                    SGPP::base::DataVector& value) override;
  virtual void clone(std::unique_ptr<VectorFunction>& clone) const override;
};

class G8ConstraintGradient : public VectorFunctionGradient {
 public:
  G8ConstraintGradient();
  virtual ~G8ConstraintGradient() override;
  virtual void eval(const SGPP::base::DataVector& x,
                    SGPP::base::DataVector& value,
                    SGPP::base::DataMatrix& gradient) override;
  virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override;
};

#endif /* OBJECTIVE_FUNCTIONS_HPP */
