// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <fstream>
#include <memory>
#include <string>
#include <algorithm>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLDevice.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingModOCLMaskMultiPlatform {

template <typename T>
class SourceBuilderMult : public base::KernelSourceBuilderBase<T> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;
  size_t dataBlockSize;
  size_t maxDimUnroll;
  size_t prefetchSize;

  std::string getData(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      output << "data_" << dataBlockingIndex << "[" << dim << "]";
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      output << "data_" << dataBlockingIndex << "_" << dim;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim << ") + "
             << dataBlockingIndex << "]";
    } else {
      throw new base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
  }

  std::string getData(size_t dim, size_t dataBlockingIndex) {
    std::stringstream dimStringStream;
    dimStringStream << dim;
    std::string dimString = dimStringStream.str();
    return this->getData(dimString, dataBlockingIndex);
  }

  std::string unrolledBasisFunctionEvalulation1D(size_t dims, size_t startDim, size_t endDim,
                                                 std::string unrollVariable) {
    std::stringstream output;

    for (size_t d = startDim; d < endDim; d++) {
      std::stringstream dimElement;
      dimElement << "(";
      if (!unrollVariable.compare("") == 0) {
        dimElement << unrollVariable << " + ";
      }
      dimElement << d;
      dimElement << ")";
      std::string pointerAccess = dimElement.str();

      std::string dString;
      if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
        std::stringstream stream;
        stream << (d);
        dString = stream.str();
      } else {
        dString = pointerAccess;
      }

      std::stringstream levelAccessStream;
      std::stringstream indexAccessStream;
      std::stringstream maskAccessStream;
      std::stringstream offsetAccessStream;
      if (useLocalMemory) {
        levelAccessStream << "locLevel[dimLevelIndex]";
        indexAccessStream << "locIndex[dimLevelIndex]";
        maskAccessStream << "locMask[dimLevelIndex]";
        offsetAccessStream << "locOffset[dimLevelIndex]";
      } else {
        levelAccessStream << "ptrLevel[dimLevelIndex]";
        indexAccessStream << "ptrIndex[dimLevelIndex]";
        maskAccessStream << "ptrMask[dimLevelIndex]";
        offsetAccessStream << "ptrOffset[dimLevelIndex]";
      }
      std::string levelAccess = levelAccessStream.str();
      std::string indexAccess = indexAccessStream.str();
      std::string maskAccess = maskAccessStream.str();
      std::string offsetAccess = offsetAccessStream.str();

      //      output << this->indent[2] << "dimLevelIndex = "
      //             << "(k * " << dims << ") + " << pointerAccess << ";" << std::endl;
      //      for (size_t i = 0; i < dataBlockSize; i++) {
      //        output << this->indent[2] << "curSupport_" << i << " *= fmax(1.0" <<
      //        this->constSuffix()
      //               << " - fabs((";
      //        output << levelAccess << " * " << getData(dString, i) << ") - " << indexAccess <<
      //        "), 0.0"
      //               << this->constSuffix() << ");" << std::endl;
      //      }

      for (size_t i = 0; i < dataBlockSize; i++) {
        output << this->indent[2] << "dimLevelIndex = "
               << "(k * " << dims << ") + " << pointerAccess << ";" << std::endl;
        output << this->indent[2] << "eval = " << levelAccess << " * " << getData(dString, i) << ";"
               << std::endl;
        output << this->indent[2] << "index_calc = eval - " << indexAccess << ";" << std::endl;
        output << this->indent[2] << "abs = as_" << this->floatType() << "(as_" << this->intType()
               << "(index_calc) | as_" << this->intType() << "(" << maskAccess << "));"
               << std::endl;
        output << this->indent[2] << "last = " << offsetAccess << " + abs;" << std::endl;
        output << this->indent[2] << "localSupport = fmax(last, 0.0" << this->constSuffix() << ");"
               << std::endl;
        output << this->indent[2] << "curSupport_" << i << " *= localSupport;" << std::endl
               << std::endl;
      }
    }
    return output.str();
  }

  std::string unrolledBasisFunctionEvalulation() {
    std::stringstream output;

    if (dims > maxDimUnroll) {
      output << this->indent[1] << "for (int unrollDim = 0; unrollDim < "
             << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll << ") {"
             << std::endl;
      output << this->unrolledBasisFunctionEvalulation1D(dims, 0, std::min(maxDimUnroll, dims),
                                                         "unrollDim");
      output << this->indent[1] << "}" << std::endl;

      if (dims % maxDimUnroll != 0) {
        output << this->unrolledBasisFunctionEvalulation1D(
            dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "");
      }
    } else {
      output << this->unrolledBasisFunctionEvalulation1D(dims, 0, dims, "");
    }

    return output.str();
  }

 public:
  SourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration,
                    size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt();
    maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
    prefetchSize = kernelConfiguration["KERNEL_PREFETCH_SIZE"].getUInt();
  }

  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("streamingModOCLMaskMP_mult.cl");
    }

    std::stringstream sourceStream;

    sourceStream << "// platform: " << device->platformName << " device: " << device->deviceName
                 << std::endl
                 << std::endl;

    if (std::is_same<T, double>::value) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl
                   << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl;
    sourceStream << "void multOCLMask(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrMask," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrOffset,"
                 << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrAlpha," << std::endl;
    sourceStream << "           __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    sourceStream << "           int resultSize," << std::endl;
    sourceStream << "           int start_grid," << std::endl;
    sourceStream << "           int end_grid) " << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << this->indent[0] << "int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int localIdx = get_local_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int groupSize = get_local_size(0);" << std::endl;
    sourceStream << this->indent[0] << "int globalSize = get_global_size(0);" << std::endl;
    sourceStream << std::endl;

    if (useLocalMemory) {
      sourceStream << "   __local " << this->floatType() << " locLevel[" << dims * prefetchSize
                   << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locIndex[" << dims * prefetchSize
                   << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locMask[" << dims * prefetchSize
                   << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locOffset[" << dims * prefetchSize
                   << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locAlpha[" << prefetchSize << "];"
                   << std::endl;
      sourceStream << std::endl;
    }

    sourceStream << "   " << this->floatType() << " eval, index_calc, abs, last, localSupport;"
                 << std::endl
                 << std::endl;

    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[0] << this->floatType() << " myResult_" << i << " = 0.0;"
                   << std::endl;
    }
    sourceStream << std::endl;

    // caching data in register array, this also requires loading the data into
    // the registers (in contrast using pointers to data directly)
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[0] << this->floatType() << " data_" << i << "[" << dims << "];"
                     << std::endl;
      }
      sourceStream << std::endl;
      for (size_t i = 0; i < dataBlockSize; i++) {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[0] << getData(d, i) << " = ptrData[(resultSize * " << d
                       << ") + (globalSize * " << i << ") + globalIdx"
                       << "];" << std::endl;
        }
        sourceStream << std::endl;
      }
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      for (size_t i = 0; i < dataBlockSize; i++) {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[0] << this->floatType() << " " << getData(d, i)
                       << " = ptrData[(resultSize * " << d << ") + (globalSize * " << i
                       << ") + globalIdx"
                       << "];" << std::endl;
        }
        sourceStream << std::endl;
      }
    }

    sourceStream << this->indent[0] << "int dimLevelIndex;" << std::endl;
    sourceStream << std::endl;

    if (useLocalMemory) {
      sourceStream << this->indent[0] << "for(int j = start_grid; j < end_grid; j+=" << prefetchSize
                   << ") {" << std::endl;

      sourceStream << this->indent[1] << "if (localIdx < " << prefetchSize << ") {" << std::endl;
      for (size_t d = 0; d < dims; d++) {
        sourceStream << this->indent[1] << "locLevel[(localIdx*" << dims << ")+" << d
                     << "] = ptrLevel[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
        sourceStream << this->indent[1] << "locIndex[(localIdx*" << dims << ")+" << d
                     << "] = ptrIndex[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
        sourceStream << this->indent[1] << "locMask[(localIdx*" << dims << ")+" << d
                     << "] = ptrMask[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
        sourceStream << this->indent[1] << "locOffset[(localIdx*" << dims << ")+" << d
                     << "] = ptrOffset[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
      }
      sourceStream << this->indent[1] << "locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
      sourceStream << this->indent[0] << "}" << std::endl;

      sourceStream << this->indent[0] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << std::endl;
      sourceStream << this->indent[0] << "for(int k = 0; k < " << prefetchSize << "; k++) {"
                   << std::endl;

      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[1] << this->floatType() << " curSupport_" << i
                     << " = locAlpha[k];" << std::endl;
      }
      sourceStream << std::endl;

      sourceStream << this->unrolledBasisFunctionEvalulation();

      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[1] << "myResult_" << i << " += curSupport_" << i << ";"
                     << std::endl;
      }
      sourceStream << this->indent[0] << "}" << std::endl;
      sourceStream << std::endl;

      sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << "   }" << std::endl;
      sourceStream << std::endl;
    } else {
      sourceStream << " for(int k = start_grid; k < end_grid; k++)" << std::endl;
      sourceStream << "   {" << std::endl;

      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[1] << this->floatType() << " curSupport_" << i
                     << " = ptrAlpha[k];" << std::endl;
      }
      sourceStream << std::endl;

      sourceStream << this->unrolledBasisFunctionEvalulation();

      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[1] << "myResult_" << i << " += curSupport_" << i << ";"
                     << std::endl;
      }
      sourceStream << this->indent[0] << "}" << std::endl;
      sourceStream << std::endl;
    }

    sourceStream << std::endl;
    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[0] << "ptrResult[(globalSize * " << i
                   << ") + globalIdx] = myResult_" << i << ";" << std::endl;
    }
    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("streamingModOCLMaskMP_mult.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};

}  // namespace StreamingModOCLMaskMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
