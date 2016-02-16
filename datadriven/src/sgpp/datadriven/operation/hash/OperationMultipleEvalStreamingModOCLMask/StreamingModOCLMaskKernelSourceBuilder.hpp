/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>
#include <memory>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

//#include "StreamingOCLParameters.hpp"

namespace SGPP {
namespace datadriven {

class StreamingModOCLMaskKernelSourceBuilder {
private:
    std::shared_ptr<base::OCLOperationConfiguration> parameters;

    size_t dims;

    size_t localWorkgroupSize;
    bool useLocalMemory;
//    size_t maxDimUnroll;

    std::string indent;

    std::string indent2;

    std::string indent3;

    std::string indent4;

    std::string reuseSource(std::string fileName);

    void writeSource(std::string fileName, std::string source);

    std::string asString();

    std::string constSuffix();

    std::string intAsString();

    void recalculateLevelIndexMask();
public:
    StreamingModOCLMaskKernelSourceBuilder(std::shared_ptr<base::OCLOperationConfiguration> parameters, size_t dims);

    std::string generateSourceMult();

    std::string generateSourceMultTrans();
};

}
}

