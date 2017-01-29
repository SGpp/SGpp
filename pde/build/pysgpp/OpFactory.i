// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "pde/src/sgpp/pde/operation/PdeOpFactory.hpp"

%newobject sgpp::op_factory::createOperationLaplace(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLaplace(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationLTwoDotProduct(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLTwoDotExplicit(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLTwoDotExplicit(
    sgpp::base::DataMatrix* m, sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLaplaceEnhanced(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLaplaceEnhanced(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
