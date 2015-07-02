/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>

#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace datadriven {

class StreamingOCLMPKernelSourceBuilder {
private:
    base::OCLConfigurationParameters parameters;
    size_t dims;

    size_t localWorkgroupSize;
    bool useLocalMemory;
    size_t dataBlockSize;
    uint64_t maxDimUnroll;

    std::string indent;

    std::string indent2;

    std::string indent3;

    std::string indent4;

    std::string asString();

    std::string constSuffix();

    std::string intAsString();

    std::string getData(std::string dim, size_t dataBlockingIndex);

    std::string getData(size_t dim, size_t dataBlockingIndex);
public:
    StreamingOCLMPKernelSourceBuilder(base::OCLConfigurationParameters parameters, size_t dims);

    std::string generateSourceMult();

    std::string generateSourceMultTrans();

    std::string reuseSource(std::string fileName);

    void writeSource(std::string fileName, std::string source);

    std::string unrolledBasisFunctionEvalulation(size_t dims, size_t startDim, size_t endDim,
            std::string unrollVariable);

}
;

}
}

