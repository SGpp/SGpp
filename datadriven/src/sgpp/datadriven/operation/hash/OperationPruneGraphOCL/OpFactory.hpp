#pragma once

#include <sgpp/globaldef.hpp>

#include "OperationPruneGraphOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationPruneGraphOCL* pruneNearestNeighborGraphConfigured(base::Grid& grid, size_t dimensions, base::DataVector &alpha,
																										 base::DataVector &data, double treshold, size_t k, std::string opencl_conf);
}
}
