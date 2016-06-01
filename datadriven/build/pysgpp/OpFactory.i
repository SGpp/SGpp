// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
    sgpp::datadriven::KernelDensityEstimator& kde);
%newobject sgpp::op_factory::createOperationInverseRosenblattTransformationKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
%newobject sgpp::op_factory::createOperationDensityMarginalizeKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
%newobject sgpp::op_factory::createOperationDensityConditionalKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
%newobject sgpp::op_factory::createOperationMultipleEval(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration);

%{
sgpp::datadriven::OperationTest*
createOperationTest(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationTest(grid).release();
}

sgpp::base::OperationMatrix*
createOperationRegularizationDiagonal(
        sgpp::base::Grid& grid, int mode, double k) {
    return sgpp::op_factory::createOperationRegularizationDiagonal(grid, mode, k).release();
}

sgpp::datadriven::OperationDensityMarginalize*
createOperationDensityMarginalize(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationDensityMarginalize(grid).release();
}

sgpp::datadriven::OperationDensityMargTo1D*
createOperationDensityMargTo1D(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationDensityMargTo1D(grid).release();
}

sgpp::datadriven::OperationDensitySampling1D*
createOperationDensitySampling1D(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationDensitySampling1D(grid).release();
}

sgpp::datadriven::OperationDensitySampling*
createOperationDensitySampling(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationDensitySampling(grid).release();
}

sgpp::datadriven::OperationDensityRejectionSampling*
createOperationDensityRejectionSampling(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationDensityRejectionSampling(grid).release();
}

sgpp::datadriven::OperationDensityConditional*
createOperationDensityConditional(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationDensityConditional(grid).release();
}

sgpp::datadriven::OperationRosenblattTransformation*
createOperationRosenblattTransformation(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationRosenblattTransformation(grid).release();
}

sgpp::datadriven::OperationTransformation1D*
createOperationRosenblattTransformation1D(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationRosenblattTransformation1D(grid).release();
}

sgpp::datadriven::OperationInverseRosenblattTransformation*
createOperationInverseRosenblattTransformation(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationInverseRosenblattTransformation(grid).release();
}

sgpp::datadriven::OperationTransformation1D*
createOperationInverseRosenblattTransformation1D(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationInverseRosenblattTransformation1D(grid).release();
}

sgpp::datadriven::OperationRosenblattTransformationKDE*
createOperationRosenblattTransformationKDE(sgpp::datadriven::KernelDensityEstimator& kde) {
    return sgpp::op_factory::createOperationRosenblattTransformationKDE(kde).release();
}

sgpp::datadriven::OperationInverseRosenblattTransformationKDE*
createOperationInverseRosenblattTransformationKDE(sgpp::datadriven::KernelDensityEstimator& kde) {
    return sgpp::op_factory::createOperationInverseRosenblattTransformationKDE(kde).release();
}

sgpp::datadriven::OperationDensityMarginalizeKDE*
createOperationDensityMarginalizeKDE(sgpp::datadriven::KernelDensityEstimator& kde) {
    return sgpp::op_factory::createOperationDensityMarginalizeKDE(kde).release();
}

sgpp::datadriven::OperationDensityConditionalKDE*
createOperationDensityConditionalKDE(sgpp::datadriven::KernelDensityEstimator& kde) {
    return sgpp::op_factory::createOperationDensityConditionalKDE(kde).release();
}

sgpp::base::OperationMultipleEval*
createOperationMultipleEval(
        sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset,
        sgpp::datadriven::OperationMultipleEvalConfiguration& configuration) {
    return sgpp::op_factory::createOperationMultipleEval(grid, dataset, configuration).release();
}
%}

sgpp::datadriven::OperationTest* createOperationTest(
    sgpp::base::Grid& grid);
sgpp::base::OperationMatrix* createOperationRegularizationDiagonal(
    sgpp::base::Grid& grid, int mode, double k);
sgpp::datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationDensitySampling* createOperationDensitySampling(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationDensityConditional* createOperationDensityConditional(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationInverseRosenblattTransformation*
createOperationInverseRosenblattTransformation(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationTransformation1D* createOperationInverseRosenblattTransformation1D(
    sgpp::base::Grid& grid);
sgpp::datadriven::OperationRosenblattTransformationKDE* createOperationRosenblattTransformationKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
sgpp::datadriven::OperationInverseRosenblattTransformationKDE*
createOperationInverseRosenblattTransformationKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
sgpp::datadriven::OperationDensityMarginalizeKDE* createOperationDensityMarginalizeKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
sgpp::datadriven::OperationDensityConditionalKDE* createOperationDensityConditionalKDE(
    sgpp::datadriven::KernelDensityEstimator& kde);
sgpp::base::OperationMultipleEval* createOperationMultipleEval(
    sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration);
