// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%newobject sgpp::op_factory::createOperationMultipleHierarchisation(
    sgpp::base::Grid& grid);

%{
sgpp::optimization::OperationMultipleHierarchisation*
createOperationMultipleHierarchisation(sgpp::base::Grid& grid) {
    return sgpp::op_factory::createOperationMultipleHierarchisation(grid).release();
}
%}

sgpp::optimization::OperationMultipleHierarchisation*
createOperationMultipleHierarchisation(sgpp::base::Grid& grid);
