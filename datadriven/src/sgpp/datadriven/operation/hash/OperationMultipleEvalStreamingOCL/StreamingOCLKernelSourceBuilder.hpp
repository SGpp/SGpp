/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace SGPP {
  namespace datadriven {

    class StreamingOCLKernelSourceBuilder {
      private:
        std::shared_ptr<base::OCLOperationConfiguration> parameters;
        size_t dims;

        size_t localWorkgroupSize;
        bool useLocalMemory;
        size_t dataBlockSize;
        size_t transGridBlockSize;
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

        std::string getDataTrans(std::string dim, size_t dataBlockingIndex);

        std::string getLevelTrans(std::string dim, size_t gridBlockingIndex);

        std::string getIndexTrans(std::string dim, size_t gridBlockingIndex);
      public:
        StreamingOCLKernelSourceBuilder(std::shared_ptr<base::OCLOperationConfiguration> parameters, size_t dims);

        std::string generateSourceMult();

        std::string generateSourceMultTrans();

        std::string reuseSource(std::string fileName);

        void writeSource(std::string fileName, std::string source);

        std::string unrolledBasisFunctionEvalulation(size_t dims, size_t startDim, size_t endDim,
            std::string unrollVariable);

        std::string unrolledBasisFunctionEvalulationTrans(size_t dims, size_t startDim, size_t endDim,
            std::string unrollVariable);

    }
    ;

  }
}

