// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject sgpp::op_factory::createOperationGamma(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationGammaLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationLB(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLE(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLD(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationLF(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationDelta(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationDeltaLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationHestonBLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationHestonCLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationHestonDLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationHestonELog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationHestonFLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationHestonGLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
%newobject sgpp::op_factory::createOperationHestonHLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationHestonKLog(
    sgpp::base::Grid& grid, double***** coef);
%newobject sgpp::op_factory::createOperationHestonX(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationHestonY(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationHestonW(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
%newobject sgpp::op_factory::createOperationHestonZ(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);

%{
sgpp::base::OperationMatrix* createOperationGamma(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationGamma(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationGammaLog(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationGammaLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationLB(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLB(grid).release();
}

sgpp::base::OperationMatrix* createOperationLE(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLE(grid).release();
}

sgpp::base::OperationMatrix* createOperationLD(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLD(grid).release();
}

sgpp::base::OperationMatrix* createOperationLF(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationLF(grid).release();
}

sgpp::base::OperationMatrix* createOperationDelta(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationDelta(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationDeltaLog(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationDeltaLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonBLog(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationHestonBLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonCLog(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationHestonCLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonDLog(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationHestonDLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonELog(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationHestonELog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonFLog(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationHestonFLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonGLog(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationHestonGLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonHLog(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationHestonHLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonKLog(
        sgpp::base::Grid& grid, double***** coef) {
    return sgpp::op_factory::createOperationHestonKLog(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonX(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationHestonX(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonY(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationHestonY(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonW(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef) {
    return sgpp::op_factory::createOperationHestonW(grid, coef).release();
}

sgpp::base::OperationMatrix* createOperationHestonZ(
        sgpp::base::Grid& grid, sgpp::base::DataVector& coef) {
    return sgpp::op_factory::createOperationHestonZ(grid, coef).release();
}
%}

sgpp::base::OperationMatrix* createOperationGamma(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationGammaLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationLB(sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLE(sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLD(sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationLF(sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationDelta(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationDeltaLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationHestonBLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationHestonCLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationHestonDLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationHestonELog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationHestonFLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationHestonGLog(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
sgpp::base::OperationMatrix* createOperationHestonHLog(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationHestonKLog(
    sgpp::base::Grid& grid, double***** coef);
sgpp::base::OperationMatrix* createOperationHestonX(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationHestonY(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationHestonW(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& coef);
sgpp::base::OperationMatrix* createOperationHestonZ(
    sgpp::base::Grid& grid, sgpp::base::DataVector& coef);
