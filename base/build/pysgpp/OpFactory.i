// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject sgpp::op_factory::createOperationHierarchisation(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationQuadrature(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationFirstMoment(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationSecondMoment(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationConvert(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationIdentity(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationEval(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationMultipleEval(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationMultipleEvalNaive(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationNaiveEval(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationNaiveEvalGradient(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationNaiveEvalHessian(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationNaiveEvalPartialDerivative(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationMakePositive(sgpp::base::Grid& grid, sgpp::base::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm, sgpp::base::MakePositiveInterpolationAlgorithm interpolationAlgorithm);

%{
sgpp::base::OperationHierarchisation* createOperationHierarchisation(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationHierarchisation(grid).release();
}

sgpp::base::OperationQuadrature* createOperationQuadrature(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationQuadrature(grid).release();
}

sgpp::base::OperationFirstMoment* createOperationFirstMoment(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationFirstMoment(grid).release();
}

sgpp::base::OperationSecondMoment* createOperationSecondMoment(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationSecondMoment(grid).release();
}

sgpp::base::OperationConvert* createOperationConvert(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationConvert(grid).release();
}

sgpp::base::OperationMatrix* createOperationIdentity(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationIdentity(grid).release();
}

sgpp::base::OperationEval* createOperationEval(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationEval(grid).release();
}

sgpp::base::OperationMultipleEval* createOperationMultipleEval(sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset) {
  return sgpp::op_factory::createOperationMultipleEval(grid, dataset).release();
}

sgpp::base::OperationMultipleEval* createOperationMultipleEvalNaive(sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset) {
  return sgpp::op_factory::createOperationMultipleEvalNaive(grid, dataset).release();
}

sgpp::base::OperationNaiveEval* createOperationNaiveEval(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationNaiveEval(grid).release();
}

sgpp::base::OperationNaiveEvalGradient* createOperationNaiveEvalGradient(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationNaiveEvalGradient(grid).release();
}

sgpp::base::OperationNaiveEvalHessian* createOperationNaiveEvalHessian(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationNaiveEvalHessian(grid).release();
}

sgpp::base::OperationNaiveEvalPartialDerivative* createOperationNaiveEvalPartialDerivative(sgpp::base::Grid& grid) {
  return sgpp::op_factory::createOperationNaiveEvalPartialDerivative(grid).release();
}

sgpp::base::OperationMakePositive* createOperationMakePositive(sgpp::base::Grid& grid, sgpp::base::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm, sgpp::base::MakePositiveInterpolationAlgorithm interpolationAlgorithm) {
   return sgpp::op_factory::createOperationMakePositive(grid, candidateSearchAlgorithm, interpolationAlgorithm).release();
}
 %}

sgpp::base::OperationHierarchisation* createOperationHierarchisation(sgpp::base::Grid& grid);
sgpp::base::OperationQuadrature* createOperationQuadrature(sgpp::base::Grid& grid);
sgpp::base::OperationFirstMoment* createOperationFirstMoment(sgpp::base::Grid& grid);
sgpp::base::OperationSecondMoment* createOperationSecondMoment(sgpp::base::Grid& grid);
sgpp::base::OperationConvert* createOperationConvert(sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationIdentity(sgpp::base::Grid& grid);
sgpp::base::OperationEval* createOperationEval(sgpp::base::Grid& grid);
sgpp::base::OperationMultipleEval* createOperationMultipleEval(sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset);
sgpp::base::OperationMultipleEval* createOperationMultipleEvalNaive(sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset);
sgpp::base::OperationNaiveEval* createOperationNaiveEval(sgpp::base::Grid& grid);
sgpp::base::OperationNaiveEvalGradient* createOperationNaiveEvalGradient(sgpp::base::Grid& grid);
sgpp::base::OperationNaiveEvalHessian* createOperationNaiveEvalHessian(sgpp::base::Grid& grid);
sgpp::base::OperationNaiveEvalPartialDerivative* createOperationNaiveEvalPartialDerivative(sgpp::base::Grid& grid);
sgpp::base::OperationMakePositive* createOperationMakePositive(sgpp::base::Grid& grid, sgpp::base::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm, sgpp::base::MakePositiveInterpolationAlgorithm interpolationAlgorithm);
