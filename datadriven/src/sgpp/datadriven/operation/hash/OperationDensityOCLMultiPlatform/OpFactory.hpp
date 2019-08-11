// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef DENSITY_OPFACTORY_H
#define DENSITY_OPFACTORY_H


#include <sgpp/globaldef.hpp>
#include <string>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OperationDensityOCLMultiPlatform.hpp>
namespace sgpp {
namespace datadriven {


DensityOCLMultiPlatform::OperationDensity*
createDensityOCLMultiPlatformConfigured(base::Grid& grid, size_t dimension,
                                        double lambda, base::OCLOperationConfiguration *parameters,
                                        size_t platform_id, size_t device_id);
/// Generates opencl density multiplication operation given opencl device and configuration file
sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity*
createDensityOCLMultiPlatformConfigured(base::Grid& grid, size_t dimension, double lambda,
                                        std::string opencl_conf, size_t platform_id,
                                        size_t device_id);
/// Generates opencl density multiplication operation given opencl device and a serialized grid
sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity*
createDensityOCLMultiPlatformConfigured(int *gridpoints, size_t gridsize, size_t dimension,
                                        double lambda, std::string opencl_conf,
                                        size_t platform_id, size_t device_id);
DensityOCLMultiPlatform::OperationDensity*
createDensityOCLMultiPlatformConfigured(int *gridpoints, size_t gridsize, size_t dimension,
                                        double lambda, base::OCLOperationConfiguration *parameters,
                                        size_t platform_id, size_t device_id);
/// Generates opencl density multiplication operation
sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity*
createDensityOCLMultiPlatformConfigured(base::Grid& grid, size_t dimension, double lambda,
                                        std::string opencl_conf);
}  // namespace datadriven
}  // namespace sgpp

#endif /* DENSITY_OPFACTORY_H */
