#pragma once

#include <sgpp/globaldef.hpp>

#include "OperationPruneGraphOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationPruneGraphOCL* pruneNearestNeighborGraphConfigured(size_t k, size_t dimensions, std::string opencl_conf);
}
}
