// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CHECK_EQUAL_FUNCTION_HPP
#define CHECK_EQUAL_FUNCTION_HPP

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/VectorFunctionHessian.hpp>

void checkEqualFunction(sgpp::base::ScalarFunction& f, sgpp::base::ScalarFunction& g);

void checkEqualFunction(sgpp::base::ScalarFunctionGradient& f,
                        sgpp::base::ScalarFunctionGradient& g);

void checkEqualFunction(sgpp::base::ScalarFunctionHessian& f, sgpp::base::ScalarFunctionHessian& g);

void checkEqualFunction(sgpp::base::VectorFunction& f, sgpp::base::VectorFunction& g);

void checkEqualFunction(sgpp::base::VectorFunctionGradient& f,
                        sgpp::base::VectorFunctionGradient& g);

void checkEqualFunction(sgpp::base::VectorFunctionHessian& f, sgpp::base::VectorFunctionHessian& g);

#endif /* CHECK_EQUAL_FUNCTION_HPP */
