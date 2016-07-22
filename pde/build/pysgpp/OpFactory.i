// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject sgpp::op_factory::createOperationLaplace(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLaplace(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationLaplaceExplicit(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLaplaceExplicit(
    sgpp::base::DataMatrix* m, sgpp::base::Grid& grid);
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

%{
sgpp::base::OperationMatrix* createOperationLaplace(
        sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLaplace(grid).release();
}

sgpp::base::OperationMatrix* createOperationLaplace(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationLaplace(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationLaplaceExplicit(
        sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLaplaceExplicit(grid).release();
}

sgpp::base::OperationMatrix* createOperationLaplaceExplicit(
        sgpp::base::DataMatrix* m, sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLaplaceExplicit(m, grid).release();
}

sgpp::base::OperationMatrix* createOperationLTwoDotProduct(
        sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLTwoDotProduct(grid).release();
}

sgpp::base::OperationMatrix* createOperationLTwoDotExplicit(
        sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLTwoDotExplicit(grid).release();
}

sgpp::base::OperationMatrix* createOperationLTwoDotExplicit(
        sgpp::base::DataMatrix* m, sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLTwoDotExplicit(m, grid).release();
}

sgpp::base::OperationMatrix* createOperationLaplaceEnhanced(
        sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLaplaceEnhanced(grid).release();
}

sgpp::base::OperationMatrix* createOperationLaplaceEnhanced(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationLaplaceEnhanced(grid, coef).release();
}
%}

sgpp::base::OperationMatrix* createOperationLaplace(
    sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLaplace(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationLaplaceExplicit(
    sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLaplaceExplicit(
    sgpp::base::DataMatrix* m, sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLTwoDotProduct(
    sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLTwoDotExplicit(
    sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLTwoDotExplicit(
    sgpp::base::DataMatrix* m, sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLaplaceEnhanced(
    sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLaplaceEnhanced(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
