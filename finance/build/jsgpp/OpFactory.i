// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject SGPP::op_factory::createOperationGamma(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationGammaLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationLB(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLE(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLD(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationLF(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationDelta(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationDeltaLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationHestonBLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationHestonCLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationHestonDLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationHestonELog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationHestonFLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationHestonGLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
%newobject SGPP::op_factory::createOperationHestonHLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationHestonKLog(
    SGPP::base::Grid& grid, SGPP::float_t***** coef);
%newobject SGPP::op_factory::createOperationHestonX(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationHestonY(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationHestonW(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
%newobject SGPP::op_factory::createOperationHestonZ(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);

%{
SGPP::base::OperationMatrix* createOperationGamma(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationGamma(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationGammaLog(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationGammaLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationLB(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLB(grid).release();
}

SGPP::base::OperationMatrix* createOperationLE(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLE(grid).release();
}

SGPP::base::OperationMatrix* createOperationLD(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLD(grid).release();
}

SGPP::base::OperationMatrix* createOperationLF(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationLF(grid).release();
}

SGPP::base::OperationMatrix* createOperationDelta(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationDelta(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationDeltaLog(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationDeltaLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonBLog(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationHestonBLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonCLog(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationHestonCLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonDLog(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationHestonDLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonELog(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationHestonELog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonFLog(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationHestonFLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonGLog(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationHestonGLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonHLog(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationHestonHLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonKLog(
        SGPP::base::Grid& grid, SGPP::float_t***** coef) {
    return SGPP::op_factory::createOperationHestonKLog(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonX(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationHestonX(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonY(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationHestonY(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonW(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef) {
    return SGPP::op_factory::createOperationHestonW(grid, coef).release();
}

SGPP::base::OperationMatrix* createOperationHestonZ(
        SGPP::base::Grid& grid, SGPP::base::DataVector& coef) {
    return SGPP::op_factory::createOperationHestonZ(grid, coef).release();
}
%}

SGPP::base::OperationMatrix* createOperationGamma(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationGammaLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationLB(SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLE(SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLD(SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationLF(SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationDelta(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationDeltaLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationHestonBLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationHestonCLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationHestonDLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationHestonELog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationHestonFLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationHestonGLog(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
SGPP::base::OperationMatrix* createOperationHestonHLog(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationHestonKLog(
    SGPP::base::Grid& grid, SGPP::float_t***** coef);
SGPP::base::OperationMatrix* createOperationHestonX(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationHestonY(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationHestonW(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& coef);
SGPP::base::OperationMatrix* createOperationHestonZ(
    SGPP::base::Grid& grid, SGPP::base::DataVector& coef);
