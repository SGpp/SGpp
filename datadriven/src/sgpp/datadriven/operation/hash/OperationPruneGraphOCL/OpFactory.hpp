// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <string>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OperationPruneGraphOCLMultiPlatform.hpp>
namespace sgpp {
namespace datadriven {

/// Generates the graph pruning operation for a specific opencl device
sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(base::Grid& grid, size_t dimensions, base::DataVector &alpha,
                                    base::DataMatrix &data, double treshold, size_t k,
                                    std::string opencl_conf, size_t platformid, size_t deviceid);
/// Generates the graph pruning operation for a specific opencl device using a serialized grid
sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(int *gridpoints, size_t gridsize, size_t dimensions,
                                    double *alpha, base::DataMatrix &data, double treshold,
                                    size_t k, std::string opencl_conf, size_t platformid,
                                    size_t deviceid);
DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(int *gridpoints, size_t gridsize, size_t dimensions,
                                    double *alpha, base::DataMatrix &data,
                                    double treshold, size_t k,
                                    sgpp::base::OCLOperationConfiguration *parameters,
                                    size_t platformid, size_t deviceid);
/// Generates the graph pruning operation
sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(base::Grid& grid, size_t dimensions, base::DataVector &alpha,
                                    base::DataMatrix &data, double treshold, size_t k,
                                    std::string opencl_conf);
}  // namespace datadriven
}  // namespace sgpp
