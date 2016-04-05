// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CREATE_GRAPH_OPFACTOR_H
#define CREATE_GRAPH_OPFACTOR_H


#include <sgpp/globaldef.hpp>
#include <string>
#include "OperationCreateGraphOCLSingleDevice.hpp"
namespace sgpp {
namespace datadriven {

DensityOCLMultiPlatform::OperationCreateGraphOCL*
createNearestNeighborGraphConfigured(base::DataMatrix &dataset, size_t k,
                                     size_t dimensions, std::string opencl_conf,
                                     size_t platformid, size_t devicdeid);
DensityOCLMultiPlatform::OperationCreateGraphOCL*
createNearestNeighborGraphConfigured(double *dataset, size_t dataset_size, size_t k,
                                     size_t dimensions, std::string opencl_conf,
                                     size_t platformid, size_t devicdeid);
}  // namespace datadriven
}  // namespace sgpp

#endif /* CREATE_GRAPH_OPFACTOR_H */
