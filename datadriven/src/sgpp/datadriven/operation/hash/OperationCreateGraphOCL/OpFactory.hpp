// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <string>
#include "OperationCreateGraphOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

DensityOCLMultiPlatform::OperationCreateGraphOCL*
createNearestNeighborGraphConfigured(base::DataMatrix &dataset, size_t k,
                                     size_t dimensions, std::string opencl_conf);
}
}
