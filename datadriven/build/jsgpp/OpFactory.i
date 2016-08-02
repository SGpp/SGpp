// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "datadriven/src/sgpp/datadriven/DatadrivenOpFactory.hpp"

%newobject sgpp::op_factory::createOperationTest(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationRegularizationDiagonal(
    sgpp::base::Grid& grid, int mode, double k);
%newobject sgpp::op_factory::createOperationDensityMarginalize(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationDensityMargTo1D(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationDensitySampling1D(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationDensitySampling(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationDensityRejectionSampling(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationDensityConditional(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationRosenblattTransformation(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationRosenblattTransformation1D(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationInverseRosenblattTransformation(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationInverseRosenblattTransformation1D(
    sgpp::base::Grid& grid);
%newobject sgpp::op_factory::createOperationRosenblattTransformationKDE(
    sgpp::datadriven::GaussianKDE& kde);
%newobject sgpp::op_factory::createOperationInverseRosenblattTransformationKDE(
    sgpp::datadriven::GaussianKDE& kde);
%newobject sgpp::op_factory::createOperationDensityMarginalizeKDE(
    sgpp::datadriven::GaussianKDE& kde);
%newobject sgpp::op_factory::createOperationDensityConditionalKDE(
    sgpp::datadriven::GaussianKDE& kde);
%newobject sgpp::op_factory::createOperationMultipleEval(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration);
