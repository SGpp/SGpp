// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject SGPP::op_factory::createOperationMultipleHierarchisation(
    SGPP::base::Grid& grid);

%{
SGPP::optimization::OperationMultipleHierarchisation*
createOperationMultipleHierarchisation(SGPP::base::Grid& grid) {
    return SGPP::op_factory::createOperationMultipleHierarchisation(grid).release();
}
%}

SGPP::optimization::OperationMultipleHierarchisation*
createOperationMultipleHierarchisation(SGPP::base::Grid& grid);
