// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject SGPP::op_factory::createOperationHierarchisation(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationQuadrature(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationFirstMoment(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationSecondMoment(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationConvert(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationIdentity(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationEval(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationMultipleEval(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationNaiveEval(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationNaiveEvalGradient(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationNaiveEvalHessian(SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationNaiveEvalPartialDerivative(SGPP::base::Grid& grid);

%{
SGPP::base::OperationHierarchisation* createOperationHierarchisation(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationHierarchisation(grid).release();
}

SGPP::base::OperationQuadrature* createOperationQuadrature(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationQuadrature(grid).release();
}

SGPP::base::OperationFirstMoment* createOperationFirstMoment(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationFirstMoment(grid).release();
}

SGPP::base::OperationSecondMoment* createOperationSecondMoment(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationSecondMoment(grid).release();
}

SGPP::base::OperationConvert* createOperationConvert(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationConvert(grid).release();
}

SGPP::base::OperationMatrix* createOperationIdentity(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationIdentity(grid).release();
}

SGPP::base::OperationEval* createOperationEval(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationEval(grid).release();
}

SGPP::base::OperationMultipleEval* createOperationMultipleEval(SGPP::base::Grid& grid, SGPP::base::DataMatrix& dataset) {
  return SGPP::op_factory::createOperationMultipleEval(grid, dataset).release();
}

SGPP::base::OperationNaiveEval* createOperationNaiveEval(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationNaiveEval(grid).release();
}

SGPP::base::OperationNaiveEvalGradient* createOperationNaiveEvalGradient(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationNaiveEvalGradient(grid).release();
}

SGPP::base::OperationNaiveEvalHessian* createOperationNaiveEvalHessian(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationNaiveEvalHessian(grid).release();
}

SGPP::base::OperationNaiveEvalPartialDerivative* createOperationNaiveEvalPartialDerivative(SGPP::base::Grid& grid) {
  return SGPP::op_factory::createOperationNaiveEvalPartialDerivative(grid).release();
}
%}

SGPP::base::OperationHierarchisation* createOperationHierarchisation(SGPP::base::Grid& grid);
SGPP::base::OperationQuadrature* createOperationQuadrature(SGPP::base::Grid& grid);
SGPP::base::OperationFirstMoment* createOperationFirstMoment(SGPP::base::Grid& grid);
SGPP::base::OperationSecondMoment* createOperationSecondMoment(SGPP::base::Grid& grid);
SGPP::base::OperationConvert* createOperationConvert(SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationIdentity(SGPP::base::Grid& grid);
SGPP::base::OperationEval* createOperationEval(SGPP::base::Grid& grid);
SGPP::base::OperationMultipleEval* createOperationMultipleEval(SGPP::base::Grid& grid, SGPP::base::DataMatrix& dataset);
SGPP::base::OperationNaiveEval* createOperationNaiveEval(SGPP::base::Grid& grid);
SGPP::base::OperationNaiveEvalGradient* createOperationNaiveEvalGradient(SGPP::base::Grid& grid);
SGPP::base::OperationNaiveEvalHessian* createOperationNaiveEvalHessian(SGPP::base::Grid& grid);
SGPP::base::OperationNaiveEvalPartialDerivative* createOperationNaiveEvalPartialDerivative(SGPP::base::Grid& grid);
