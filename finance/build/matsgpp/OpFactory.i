// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "finance/src/sgpp/finance/operation/FinanceOpFactory.hpp"

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
