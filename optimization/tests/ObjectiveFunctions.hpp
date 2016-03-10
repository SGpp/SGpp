// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OBJECTIVE_FUNCTIONS_HPP
#define OBJECTIVE_FUNCTIONS_HPP

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionHessian.hpp>

#include <cmath>

using sgpp::optimization::ScalarFunction;
using sgpp::optimization::ScalarFunctionGradient;
using sgpp::optimization::ScalarFunctionHessian;
using sgpp::optimization::VectorFunction;
using sgpp::optimization::VectorFunctionGradient;
using sgpp::optimization::VectorFunctionHessian;

class ExampleFunction : public ScalarFunction {
 public:
  ExampleFunction();
  ~ExampleFunction() override;
  double eval(const sgpp::base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

class ExampleGradient : public ScalarFunctionGradient {
 public:
  ExampleGradient();
  ~ExampleGradient() override;
  double eval(const sgpp::base::DataVector& x,
                     sgpp::base::DataVector& gradient) override;
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};

class ExampleHessian : public ScalarFunctionHessian {
 public:
  ExampleHessian();
  ~ExampleHessian() override;
  double eval(const sgpp::base::DataVector& x,
                     sgpp::base::DataVector& gradient,
                     sgpp::base::DataMatrix& hessian) override;
  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const
  override;
};



class SphereGradient : public ScalarFunctionGradient {
 public:
  explicit SphereGradient(size_t d);
  ~SphereGradient() override;
  double eval(const sgpp::base::DataVector& x,
                     sgpp::base::DataVector& gradient) override;
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};

class SphereHessian : public ScalarFunctionHessian {
 public:
  explicit SphereHessian(size_t d);
  ~SphereHessian() override;
  double eval(const sgpp::base::DataVector& x,
                     sgpp::base::DataVector& gradient,
                     sgpp::base::DataMatrix& hessian) override;
  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const
  override;
};



class DeformedLinearPhiFunction : public VectorFunction {
 public:
  explicit DeformedLinearPhiFunction(size_t d);
  ~DeformedLinearPhiFunction() override;
  void eval(const sgpp::base::DataVector& x,
            sgpp::base::DataVector& value) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;

 protected:
  sgpp::base::DataVector eigenvalues;
};

class DeformedLinearPhiGradient : public VectorFunctionGradient {
 public:
  explicit DeformedLinearPhiGradient(size_t d);
  ~DeformedLinearPhiGradient() override;
  void eval(const sgpp::base::DataVector& x,
            sgpp::base::DataVector& value,
            sgpp::base::DataMatrix& gradient) override;
  void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override;

 protected:
  sgpp::base::DataVector eigenvalues;
};



class G3ObjectiveFunction : public ScalarFunction {
 public:
  explicit G3ObjectiveFunction(size_t d);
  ~G3ObjectiveFunction() override;
  double eval(const sgpp::base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

class G3ObjectiveGradient : public ScalarFunctionGradient {
 public:
  explicit G3ObjectiveGradient(size_t d);
  ~G3ObjectiveGradient() override;
  double eval(const sgpp::base::DataVector& x,
                     sgpp::base::DataVector& gradient) override;
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};



class G3ConstraintFunction : public VectorFunction {
 public:
  explicit G3ConstraintFunction(size_t d);
  ~G3ConstraintFunction() override;
  void eval(const sgpp::base::DataVector& x,
            sgpp::base::DataVector& value) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;
};

class G3ConstraintGradient : public VectorFunctionGradient {
 public:
  explicit G3ConstraintGradient(size_t d);
  ~G3ConstraintGradient() override;
  void eval(const sgpp::base::DataVector& x,
            sgpp::base::DataVector& value,
            sgpp::base::DataMatrix& gradient) override;
  void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override;
};



class G8ObjectiveFunction : public ScalarFunction {
 public:
  G8ObjectiveFunction();
  ~G8ObjectiveFunction() override;
  double eval(const sgpp::base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

class G8ObjectiveGradient : public ScalarFunctionGradient {
 public:
  G8ObjectiveGradient();
  ~G8ObjectiveGradient() override;
  double eval(const sgpp::base::DataVector& x,
                     sgpp::base::DataVector& gradient) override;
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override;
};



class G8ConstraintFunction : public VectorFunction {
 public:
  G8ConstraintFunction();
  ~G8ConstraintFunction() override;
  void eval(const sgpp::base::DataVector& x,
            sgpp::base::DataVector& value) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;
};

class G8ConstraintGradient : public VectorFunctionGradient {
 public:
  G8ConstraintGradient();
  ~G8ConstraintGradient() override;
  void eval(const sgpp::base::DataVector& x,
            sgpp::base::DataVector& value,
            sgpp::base::DataMatrix& gradient) override;
  void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override;
};

#endif /* OBJECTIVE_FUNCTIONS_HPP */
