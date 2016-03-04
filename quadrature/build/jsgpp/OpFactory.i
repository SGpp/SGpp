// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject sgpp::op_factory::createOperationQuadratureMCAdvanced(
    base::Grid& grid, size_t numberOfSamples, std::uint64_t seed);

%{
sgpp::quadrature::OperationQuadratureMCAdvanced*
createOperationQuadratureMCAdvanced(
        sgpp::base::Grid& grid, size_t numberOfSamples,
        std::uint64_t seed = std::mt19937_64::default_seed) {
    return sgpp::op_factory::createOperationQuadratureMCAdvanced(grid, numberOfSamples, seed).release();
}
%}

sgpp::quadrature::OperationQuadratureMCAdvanced*
createOperationQuadratureMCAdvanced(
    sgpp::base::Grid& grid, size_t numberOfSamples,
    std::uint64_t seed = std::mt19937_64::default_seed);
