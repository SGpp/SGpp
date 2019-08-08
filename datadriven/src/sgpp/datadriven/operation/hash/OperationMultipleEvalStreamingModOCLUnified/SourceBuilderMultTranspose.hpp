// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>
#include <sgpp/base/opencl/OCLDevice.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingModOCLUnified {

template <typename T>
class SourceBuilderMultTranspose : public base::KernelSourceBuilderBase<T> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;
  size_t transGridBlockSize;
  uint64_t maxDimUnroll;
  size_t transPrefetchSize;

  std::string getLevel(std::string dim, size_t gridBlockingIndex) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      output << "level_" << gridBlockingIndex << "[" << dim << "]";
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      output << "level_" << gridBlockingIndex << "_" << dim;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      output << "ptrLevel[(((globalSize * " << gridBlockingIndex << ") + globalIdx) * " << dims
             << ") + " << dim << "]";
    } else {
      throw new base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
  }

  std::string getLevel(size_t dim, size_t dataBlockingIndex) {
    return this->getLevel(std::to_string(dim), dataBlockingIndex);
  }

  std::string getIndex(std::string dim, size_t gridBlockingIndex) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      output << "index_" << gridBlockingIndex << "[" << dim << "]";
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      output << "index_" << gridBlockingIndex << "_" << dim;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      output << "ptrIndex[(((globalSize * " << gridBlockingIndex << ") + globalIdx) * " << dims
             << ") + " << dim << "]";
    } else {
      throw new base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
  }

  std::string getIndex(size_t dim, size_t dataBlockingIndex) {
    return this->getIndex(std::to_string(dim), dataBlockingIndex);
  }

  std::string getData(std::string dim) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool()) {
      output << "locData[(" << dim << " * " << transPrefetchSize << ") + k]";
    } else {
      output << "ptrData[(" << dim << " * rangeData) + k]";
    }
    return output.str();
  }

  std::string getData(size_t dim) { return this->getData(std::to_string(dim)); }

  std::string unrolledBasisFunctionEvalulation(size_t dims, size_t startDim, size_t endDim,
                                               std::string unrollVariable, size_t indentLevel) {
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

      for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
        // TODO(pfandedd): add 1d blocked basis function evaluation

        //        output << this->indent[2] << "// nothing to do on l == 1" << std::endl;
        //        output << this->indent[2] << "if (" << getLevel(dString, gridIndex) << " == 0) {"
        //               << std::endl;
        //        output << this->indent[2] << "} else {" << std::endl;
        output << this->indent[2] << "curSupport_" << gridIndex << " *= fmax(1.0"
               << this->constSuffix() << " - fabs((";
        output << getLevel(dString, gridIndex) << " * " << getData(dString) << ") - "
               << getIndex(dString, gridIndex) << "), 0.0" << this->constSuffix() << ");"
               << std::endl;

        //        output << this->indent[3] << "if (" << getIndex(dString, gridIndex) << " == 0 || "
        //               << getIndex(dString, gridIndex) << " == " << getLevel(dString, gridIndex)
        //               << ") {"
        //               << std::endl;
        //        output << this->indent[4] << "curSupport_" << gridIndex << " *= 2;" << std::endl;
        //        output << this->indent[3] << "}" << std::endl;
        //        output << this->indent[2] << "}" << std::endl;

        //                output << this->indent[2] << "} else {" << std::endl;
        //                output << this->indent[3] << "curSupport_" << i << " *= fmax(1.0" <<
        //                this->constSuffix()
        //                       << " - fabs((";
        //                output << levelAccess << " * " << getData(dString, i) << ") - " <<
        //                indexAccess <<
        //                "), 0.0"
        //                       << this->constSuffix() << ");" << std::endl;
        //                output << this->indent[3] << "if (" << indexAccess << " == 0 || " <<
        //                indexAccess
        //                       << " == " << levelAccess << ") {" << std::endl;
        //                output << this->indent[4] << "curSupport_" << i << " *= 2;" << std::endl;
        //                output << this->indent[3] << "}" << std::endl;
        //                output << this->indent[2] << "}" << std::endl;
      }
    }
    return output.str();
  }

 public:
  SourceBuilderMultTranspose(std::shared_ptr<base::OCLDevice> device,
                             json::Node &kernelConfiguration, size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    transGridBlockSize = kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();
    maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
    transPrefetchSize = kernelConfiguration["KERNEL_TRANS_PREFETCH_SIZE"].getUInt();
  }

  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("streamingModOCLUnified_multTranspose.cl");
    }

    std::stringstream sourceStream;

    size_t indentLevel = 0;

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
    sourceStream << "void multTransOCLUnified(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrSource,"
                 << std::endl;
    sourceStream << "           __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    sourceStream << "           int dataBlockStart," << std::endl;
    sourceStream << "           int dataBlockEnd)" << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << this->indent[indentLevel] << "int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << this->indent[indentLevel] << "int localIdx = get_local_id(0);" << std::endl;
    sourceStream << this->indent[indentLevel] << "int groupSize = get_local_size(0);" << std::endl;
    sourceStream << this->indent[indentLevel] << "int globalSize = get_global_size(0);"
                 << std::endl;
    sourceStream << this->indent[indentLevel] << "int rangeData = dataBlockEnd - dataBlockStart;"
                 << std::endl;

    sourceStream << std::endl;
    sourceStream << this->indent[indentLevel] << this->floatType()
                 << " eval, index_calc, abs, last, localSupport;" << std::endl
                 << std::endl;

    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
      sourceStream << this->indent[indentLevel] << this->floatType() << " myResult_" << gridIndex
                   << " = 0.0;" << std::endl;
    }
    sourceStream << std::endl;

    if (useLocalMemory) {
      sourceStream << this->indent[indentLevel] << "__local " << this->floatType() << " locData["
                   << dims * transPrefetchSize << "];" << std::endl;
      sourceStream << this->indent[indentLevel] << "__local " << this->floatType() << " locSource["
                   << transPrefetchSize << "];" << std::endl
                   << std::endl;
    }

    // create a register storage for the level and index of the grid points of
    // the work item
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
        sourceStream << this->indent[indentLevel] << this->floatType() << " level_" << gridIndex
                     << "[" << dims << "];" << std::endl;
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[indentLevel] << "level_" << gridIndex << "[" << d
                       << "] = ptrLevel[(((globalSize * " << gridIndex << ") + globalIdx) * "
                       << dims << ") + " << d << "];" << std::endl;
        }
        sourceStream << std::endl;

        sourceStream << this->indent[indentLevel] << this->floatType() << " index_" << gridIndex
                     << "[" << dims << "];" << std::endl;
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[indentLevel] << "index_" << gridIndex << "[" << d
                       << "] = ptrIndex[(((globalSize * " << gridIndex << ") + globalIdx) * "
                       << dims << ") + " << d << "];" << std::endl;
        }
        sourceStream << std::endl;
      }
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[indentLevel] << this->floatType() << " level_" << gridIndex
                       << "_" << d << " = ptrLevel[(((globalSize * " << gridIndex
                       << ") + globalIdx) * " << dims << ") + " << d << "];" << std::endl;
        }
        sourceStream << std::endl;

        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[indentLevel] << this->floatType() << " index_" << gridIndex
                       << "_" << d << " = ptrIndex[(((globalSize * " << gridIndex
                       << ") + globalIdx) * " << dims << ") + " << d << "];" << std::endl;
        }
        sourceStream << std::endl;
      }
    }

    sourceStream << std::endl;
    sourceStream << this->indent[indentLevel] << "// Iterate over all grid points" << std::endl;
    if (useLocalMemory) {
      sourceStream << this->indent[indentLevel]
                   << "for(int i = dataBlockStart; i < dataBlockEnd; i+=" << transPrefetchSize
                   << ") {" << std::endl;

      indentLevel += 1;

      sourceStream << this->indent[indentLevel] << "if (localIdx < " << transPrefetchSize << ") {"
                   << std::endl;
      indentLevel += 1;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << this->indent[indentLevel] << "locData[(" << d << "*" << transPrefetchSize
                     << ")+(localIdx)] = ptrData[(" << d << "*rangeData)+(localIdx+i)];"
                     << std::endl;
      }

      sourceStream << this->indent[indentLevel] << "locSource[localIdx] = ptrSource[i+localIdx];"
                   << std::endl;
      indentLevel -= 1;
      sourceStream << this->indent[indentLevel] << "}" << std::endl;
      sourceStream << this->indent[indentLevel] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << std::endl;
      sourceStream << this->indent[indentLevel] << "for(int k = 0; k < " << transPrefetchSize
                   << "; k++) {" << std::endl;

      indentLevel += 1;

      for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
        sourceStream << this->indent[indentLevel] << this->floatType() << " curSupport_"
                     << gridIndex << " = locSource[k];" << std::endl;
      }
    } else {
      sourceStream << this->indent[indentLevel]
                   << "for(int k = dataBlockStart; k < dataBlockEnd; k++) {" << std::endl;

      indentLevel += 1;

      for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
        sourceStream << this->indent[indentLevel] << this->floatType() << " curSupport_"
                     << gridIndex << " = ptrSource[k];" << std::endl;
      }
    }

    sourceStream << std::endl;

    if (dims > maxDimUnroll) {
      sourceStream << this->indent[indentLevel] << "for (int unrollDim = 0; unrollDim < "
                   << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll
                   << ") {" << std::endl;

      indentLevel += 1;

      sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, std::min(maxDimUnroll, dims),
                                                             "unrollDim", indentLevel);

      indentLevel -= 1;

      sourceStream << this->indent[indentLevel] << "}" << std::endl;

      if (dims % maxDimUnroll != 0) {
        sourceStream << this->unrolledBasisFunctionEvalulation(
            dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "", indentLevel);
      }

    } else {
      sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "", indentLevel);
    }

    sourceStream << std::endl;

    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
      sourceStream << this->indent[indentLevel] << "myResult_" << gridIndex << " += curSupport_"
                   << gridIndex << ";" << std::endl;
    }

    indentLevel -= 1;

    sourceStream << this->indent[indentLevel] << "}" << std::endl
                 << std::endl;

    if (useLocalMemory) {
      sourceStream << this->indent[indentLevel] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;

      indentLevel -= 1;

      sourceStream << this->indent[indentLevel] << "}" << std::endl;
    }

    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
      sourceStream << this->indent[indentLevel] << "ptrResult[(globalSize * " << gridIndex
                   << ") + globalIdx] = myResult_" << gridIndex << ";" << std::endl;
    }

    indentLevel -= 1;

    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("streamingModOCLUnified_multTranspose.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};

}  // namespace StreamingModOCLUnified
}  // namespace datadriven
}  // namespace sgpp
