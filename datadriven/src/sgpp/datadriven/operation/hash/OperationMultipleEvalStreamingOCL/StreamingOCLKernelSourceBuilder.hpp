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

class StreamingOCLKernelSourceBuilder {
private:
    base::OCLConfigurationParameters parameters;

    std::string asString();
    std::string constSuffix();
    std::string intAsString();
public:
    StreamingOCLKernelSourceBuilder(base::OCLConfigurationParameters parameters);

    std::string generateSourceMult(size_t dims);

    std::string generateSourceMultTrans(size_t dims);

    std::string reuseSource(std::string fileName);

    void writeSource(std::string fileName, std::string source);

}
;

}
}

