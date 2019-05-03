// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "quadrature/src/sgpp/quadrature/QuadratureOpFactory.hpp"

%newobject sgpp::op_factory::createOperationQuadratureMCAdvanced(
    base::Grid& grid, size_t numberOfSamples, std::uint64_t seed);
