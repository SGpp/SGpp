// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "base/src/sgpp/base/operation/BaseOpFactory.hpp"

%newobject sgpp::op_factory::createOperationHierarchisation(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationArbitraryBoundaryHierarchisation(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationQuadrature(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationFirstMoment(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationSecondMoment(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationConvert(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationIdentity(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationEval(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationMultipleEval(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationEvalNaive(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationEvalGradientNaive(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationEvalHessianNaive(sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationEvalPartialDerivativeNaive(sgpp::base::Grid& grid);
