// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject SGPP::op_factory::createOperationLaplace(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLaplace(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationLTwoDotProduct(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLTwoDotExplicit(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLTwoDotExplicit(
    SGPP::base::DataMatrix* m, SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLaplaceEnhanced(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLaplaceEnhanced(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);

%{
SGPP::base::OperationMatrix* createOperationLaplace(
        SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLaplace(grid).release();
}

SGPP::base::OperationMatrix* createOperationLaplace(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationLaplace(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationLTwoDotProduct(
        SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLTwoDotProduct(grid).release();
}

SGPP::base::OperationMatrix* createOperationLTwoDotExplicit(
        SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLTwoDotExplicit(grid).release();
}

SGPP::base::OperationMatrix* createOperationLTwoDotExplicit(
        SGPP::base::DataMatrix* m, SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLTwoDotExplicit(m, grid).release();
}

SGPP::base::OperationMatrix* createOperationLaplaceEnhanced(
        SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLaplaceEnhanced(grid).release();
}

SGPP::base::OperationMatrix* createOperationLaplaceEnhanced(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationLaplaceEnhanced(grid, coef).release();
}
%}

SGPP::base::OperationMatrix* createOperationLaplace(
    SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLaplace(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationLTwoDotProduct(
    SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLTwoDotExplicit(
    SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLTwoDotExplicit(
    SGPP::base::DataMatrix* m, SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLaplaceEnhanced(
    SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLaplaceEnhanced(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
