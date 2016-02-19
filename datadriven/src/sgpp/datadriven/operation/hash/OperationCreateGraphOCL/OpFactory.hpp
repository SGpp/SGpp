#pragma once

#include <sgpp/globaldef.hpp>

#include "OperationCreateGraphOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL* createNearestNeighborGraphConfigured(base::DataMatrix &dataset, size_t k, size_t dimensions, std::string opencl_conf);
}
}
