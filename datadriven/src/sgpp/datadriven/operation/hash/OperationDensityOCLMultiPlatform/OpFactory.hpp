// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#pragma once

#include <sgpp/globaldef.hpp>
#include <string>
#include "OperationDensityOCLMultiPlatform.hpp"
namespace sgpp {
namespace datadriven {

sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL*
createDensityOCLMultiPlatformConfigured(base::Grid& grid, size_t dimension, double lambda,
                                        std::string opencl_conf);
sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL*
createDensityOCLMultiPlatformConfigured(int *gridpoints, size_t gridsize, size_t dimension, double lambda,
                                        std::string opencl_conf);
}  // namespace datadriven
}  // namespace sgpp
