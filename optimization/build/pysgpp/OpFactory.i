// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "optimization/src/sgpp/optimization/operation/OptimizationOpFactory.hpp"

%newobject sgpp::op_factory::createOperationMultipleHierarchisation(
    sgpp::base::Grid& grid);
