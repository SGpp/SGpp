// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject SGPP::op_factory::createOperationTest(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationRegularizationDiagonal(
    SGPP::base::Grid& grid, int mode, float_t k);
%newobject SGPP::op_factory::createOperationDensityMarginalize(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationDensityMargTo1D(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationDensitySampling1D(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationDensitySampling(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationDensityRejectionSampling(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationDensityConditional(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationRosenblattTransformation(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationRosenblattTransformation1D(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationInverseRosenblattTransformation(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationInverseRosenblattTransformation1D(
    SGPP::base::Grid& grid);
%newobject SGPP::op_factory::createOperationRosenblattTransformationKDE(
    SGPP::datadriven::GaussianKDE& kde);
%newobject SGPP::op_factory::createOperationInverseRosenblattTransformationKDE(
    SGPP::datadriven::GaussianKDE& kde);
%newobject SGPP::op_factory::createOperationDensityMarginalizeKDE(
    SGPP::datadriven::GaussianKDE& kde);
%newobject SGPP::op_factory::createOperationDensityConditionalKDE(
    SGPP::datadriven::GaussianKDE& kde);
%newobject SGPP::op_factory::createOperationMultipleEval(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& dataset,
    SGPP::datadriven::OperationMultipleEvalConfiguration& configuration);

%{
SGPP::datadriven::OperationTest*
createOperationTest(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationTest(grid).release();
}

SGPP::base::OperationMatrix*
createOperationRegularizationDiagonal(
        SGPP::base::Grid& grid, int mode, float_t k) {
    return SGPP::op_factory::createOperationRegularizationDiagonal(grid, mode, k).release();
}

SGPP::datadriven::OperationDensityMarginalize*
createOperationDensityMarginalize(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationDensityMarginalize(grid).release();
}

SGPP::datadriven::OperationDensityMargTo1D*
createOperationDensityMargTo1D(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationDensityMargTo1D(grid).release();
}

SGPP::datadriven::OperationDensitySampling1D*
createOperationDensitySampling1D(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationDensitySampling1D(grid).release();
}

SGPP::datadriven::OperationDensitySampling*
createOperationDensitySampling(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationDensitySampling(grid).release();
}

SGPP::datadriven::OperationDensityRejectionSampling*
createOperationDensityRejectionSampling(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationDensityRejectionSampling(grid).release();
}

SGPP::datadriven::OperationDensityConditional*
createOperationDensityConditional(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationDensityConditional(grid).release();
}

SGPP::datadriven::OperationRosenblattTransformation*
createOperationRosenblattTransformation(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationRosenblattTransformation(grid).release();
}

SGPP::datadriven::OperationTransformation1D*
createOperationRosenblattTransformation1D(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationRosenblattTransformation1D(grid).release();
}

SGPP::datadriven::OperationInverseRosenblattTransformation*
createOperationInverseRosenblattTransformation(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationInverseRosenblattTransformation(grid).release();
}

SGPP::datadriven::OperationTransformation1D*
createOperationInverseRosenblattTransformation1D(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationInverseRosenblattTransformation1D(grid).release();
}

SGPP::datadriven::OperationRosenblattTransformationKDE*
createOperationRosenblattTransformationKDE(SGPP::datadriven::GaussianKDE& kde) {
    return SGPP::op_factory::createOperationRosenblattTransformationKDE(kde).release();
}

SGPP::datadriven::OperationInverseRosenblattTransformationKDE*
createOperationInverseRosenblattTransformationKDE(SGPP::datadriven::GaussianKDE& kde) {
    return SGPP::op_factory::createOperationInverseRosenblattTransformationKDE(kde).release();
}

SGPP::datadriven::OperationDensityMarginalizeKDE*
createOperationDensityMarginalizeKDE(SGPP::datadriven::GaussianKDE& kde) {
    return SGPP::op_factory::createOperationDensityMarginalizeKDE(kde).release();
}

SGPP::datadriven::OperationDensityConditionalKDE*
createOperationDensityConditionalKDE(SGPP::datadriven::GaussianKDE& kde) {
    return SGPP::op_factory::createOperationDensityConditionalKDE(kde).release();
}

SGPP::base::OperationMultipleEval*
createOperationMultipleEval(
        SGPP::base::Grid& grid, SGPP::base::DataMatrix& dataset,
        SGPP::datadriven::OperationMultipleEvalConfiguration& configuration) {
    return SGPP::op_factory::createOperationMultipleEval(grid, dataset, configuration).release();
}
%}

SGPP::datadriven::OperationTest* createOperationTest(
    SGPP::base::Grid& grid);
SGPP::base::OperationMatrix* createOperationRegularizationDiagonal(
    SGPP::base::Grid& grid, int mode, float_t k);
SGPP::datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationDensitySampling* createOperationDensitySampling(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationDensityConditional* createOperationDensityConditional(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationInverseRosenblattTransformation*
createOperationInverseRosenblattTransformation(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationTransformation1D* createOperationInverseRosenblattTransformation1D(
    SGPP::base::Grid& grid);
SGPP::datadriven::OperationRosenblattTransformationKDE* createOperationRosenblattTransformationKDE(
    SGPP::datadriven::GaussianKDE& kde);
SGPP::datadriven::OperationInverseRosenblattTransformationKDE*
createOperationInverseRosenblattTransformationKDE(
    SGPP::datadriven::GaussianKDE& kde);
SGPP::datadriven::OperationDensityMarginalizeKDE* createOperationDensityMarginalizeKDE(
    SGPP::datadriven::GaussianKDE& kde);
SGPP::datadriven::OperationDensityConditionalKDE* createOperationDensityConditionalKDE(
    SGPP::datadriven::GaussianKDE& kde);
SGPP::base::OperationMultipleEval* createOperationMultipleEval(
    SGPP::base::Grid& grid, SGPP::base::DataMatrix& dataset,
    SGPP::datadriven::OperationMultipleEvalConfiguration& configuration);
