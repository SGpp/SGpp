// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CHECK_EQUAL_FUNCTION_HPP
#define CHECK_EQUAL_FUNCTION_HPP

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionHessian.hpp>

void checkEqualFunction(sgpp::optimization::ScalarFunction& f,
                        sgpp::optimization::ScalarFunction& g);

void checkEqualFunction(sgpp::optimization::ScalarFunctionGradient& f,
                        sgpp::optimization::ScalarFunctionGradient& g);

void checkEqualFunction(sgpp::optimization::ScalarFunctionHessian& f,
                        sgpp::optimization::ScalarFunctionHessian& g);

void checkEqualFunction(sgpp::optimization::VectorFunction& f,
                        sgpp::optimization::VectorFunction& g);

void checkEqualFunction(sgpp::optimization::VectorFunctionGradient& f,
                        sgpp::optimization::VectorFunctionGradient& g);

void checkEqualFunction(sgpp::optimization::VectorFunctionHessian& f,
                        sgpp::optimization::VectorFunctionHessian& g);

#endif /* CHECK_EQUAL_FUNCTION_HPP */
