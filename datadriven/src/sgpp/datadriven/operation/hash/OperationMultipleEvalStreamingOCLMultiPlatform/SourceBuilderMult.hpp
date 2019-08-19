// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <fstream>
#include <string>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

/**
 * Class for creating compute kernels for the MultiEval operation \f$v:= B^T \alpha\f$.
 * This class uses the parameters of the provided configuration to create OpenCL source code.
 * The generated code is specific to a single device, there each KernelMult requires its own
 * SourceBuilderMult.
 * For the code generation, a code fragment concatenation approach is used.
 *
 * @see StreamingOCLMultiPlatform::Configuration
 * @see StreamingOCLMultiPlatform::KernelMult
 */
template <typename real_type>
class SourceBuilderMult : public base::KernelSourceBuilderBase<real_type> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;
  size_t dataBlockSize;
  uint64_t maxDimUnroll;
  size_t prefetchSize;

  /**
   * Generates code fragments that models a dataset access.
   *
   * @param dim dimension to be accessed
   * @param dataBlockingIndex data blocking information to be taken into account
   */
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
      throw base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
  }

  /**
   * Generates code fragments that models a dataset access.
   *
   * @param dim dimension to be accessed
   * @param dataBlockingIndex data blocking information to be taken into account
   */
  std::string getData(size_t dim, size_t dataBlockingIndex) {
    std::stringstream dimStringStream;
    dimStringStream << dim;
    std::string dimString = dimStringStream.str();
    return this->getData(dimString, dataBlockingIndex);
  }

  /**
   * Generates code fragments for an basis function evaluation.
   * Creates only dimension after another without a loop.
   *
   * @param dims Number of dimensions to be implemented
   * @param startDim First dimension of the unroll, useful for implementing loop unrolling
   * @param endDim Last dimension of the unrolled implementation
   * @param unrollVariable Variable used for partial unrolling, can be an empty string if everything
   * is unrolled
   */
  std::string unrolledBasisFunctionEvalulation(size_t dims, size_t startDim, size_t endDim,
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
      if (useLocalMemory) {
        levelAccessStream << "locLevel[dimLevelIndex]";
        indexAccessStream << "locIndex[dimLevelIndex]";
      } else {
        levelAccessStream << "ptrLevel[dimLevelIndex]";
        indexAccessStream << "ptrIndex[dimLevelIndex]";
      }
      std::string levelAccess = levelAccessStream.str();
      std::string indexAccess = indexAccessStream.str();

      output << this->indent[2] << "dimLevelIndex = "
             << "(k * " << dims << ") + " << pointerAccess << ";" << std::endl;

      for (size_t i = 0; i < dataBlockSize; i++) {
        output << this->indent[2] << "curSupport_" << i << " *= fmax(1.0" << this->constSuffix()
               << " - fabs((";
        output << levelAccess << " * " << getData(dString, i) << ") - " << indexAccess << "), 0.0"
               << this->constSuffix() << ");" << std::endl;
      }
    }
    return output.str();
  }

 public:
  /**
   * Constructor
   *
   * @param device The device this code generators creates for
   * @param kernelConfiguration The configuration for this device
   * @param dims Dimension of the data mining problem
   */
  SourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration,
                    size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt();
    maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
    prefetchSize = kernelConfiguration["KERNEL_PREFETCH_SIZE"].getUInt();
  }

  /**
   * Entry point of source generator.
   * Creates compute kernel for the MultiEval operation.
   */
  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("StreamingOCLMultiPlatform_mult.cl");
    }

    std::stringstream sourceStream;

    if (std::is_same<real_type, double>::value) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl;
    sourceStream << "void multOCL(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "             __global const " << this->floatType() << "* ptrIndex,"
                 << std::endl;
    sourceStream << "             __global const " << this->floatType() << "* ptrData,"
                 << std::endl;
    sourceStream << "             __global const " << this->floatType() << "* ptrAlpha,"
                 << std::endl;
    sourceStream << "             __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    sourceStream << "             int resultSize," << std::endl;
    sourceStream << "             int start_grid," << std::endl;
    sourceStream << "             int end_grid) {" << std::endl;
    sourceStream << this->indent[0] << "int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int localIdx = get_local_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int globalSize = get_global_size(0);" << std::endl;

    sourceStream << std::endl;

    if (useLocalMemory) {
      sourceStream << this->indent[0] << "__local " << this->floatType() << " locLevel["
                   << dims * prefetchSize << "];" << std::endl;
      sourceStream << this->indent[0] << "__local " << this->floatType() << " locIndex["
                   << dims * prefetchSize << "];" << std::endl;
      sourceStream << this->indent[0] << "__local " << this->floatType() << " locAlpha["
                   << prefetchSize << "];" << std::endl;
      sourceStream << std::endl;
    }

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
      //      sourceStream << this->indent[0] << "int dataBase = "
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
      sourceStream
          << this->indent[0]
          << "for(int gridBlockStart = start_grid; gridBlockStart < end_grid; gridBlockStart += "
          << prefetchSize << ") {" << std::endl;

      sourceStream << this->indent[1] << "if (localIdx < " << prefetchSize << ") {" << std::endl;
      if (dims > maxDimUnroll) {
        sourceStream << this->indent[1] << "for (int d = 0; d < " << dims << "; d++) {"
                     << std::endl;
        sourceStream << this->indent[2] << "locLevel[d * " << prefetchSize
                     << " + localIdx] = ptrLevel[gridBlockStart * " << dims << " +  ("
                     << prefetchSize << " * d) + localIdx];" << std::endl;
        sourceStream << this->indent[1] << "}" << std::endl;

        sourceStream << this->indent[1] << "for (int d = 0; d < " << dims << "; d++) {"
                     << std::endl;
        sourceStream << this->indent[2] << "locIndex[d * " << prefetchSize
                     << " + localIdx] = ptrIndex[gridBlockStart * " << dims << " +  ("
                     << prefetchSize << " * d) + localIdx];" << std::endl;
        sourceStream << this->indent[1] << "}" << std::endl;
      } else {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[1] << "locLevel[" << d << " * " << prefetchSize
                       << " + localIdx] = ptrLevel[gridBlockStart * " << dims << " +  ("
                       << prefetchSize << " * " << d << ") + localIdx];" << std::endl;
        }
        sourceStream << this->indent[1] << std::endl;
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[1] << "locIndex[" << d << " * " << prefetchSize
                       << " + localIdx] = ptrIndex[gridBlockStart * " << dims << " +  ("
                       << prefetchSize << " * " << d << ") + localIdx];" << std::endl;
        }
      }

      sourceStream << this->indent[1] << std::endl;
      sourceStream << this->indent[1] << "locAlpha[localIdx] = ptrAlpha[gridBlockStart + localIdx];"
                   << std::endl;

      sourceStream << this->indent[1] << "}" << std::endl;
      sourceStream << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << std::endl;
      sourceStream << this->indent[1] << "for(int k = 0; k < " << prefetchSize << "; k++) {"
                   << std::endl;
      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[2] << this->floatType() << " curSupport_" << i
                     << " = locAlpha[k];" << std::endl;
      }
      sourceStream << std::endl;
    } else {
      sourceStream << this->indent[0] << "for(int k = start_grid; k < end_grid; k++) {"
                   << std::endl;
      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[1] << this->floatType() << " curSupport_" << i
                     << " = ptrAlpha[k];" << std::endl;
      }
      sourceStream << std::endl;
    }

    if (dims > maxDimUnroll) {
      sourceStream << this->indent[1] << "for (int unrollDim = 0; unrollDim < "
                   << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll
                   << ") {" << std::endl;
      sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, std::min(maxDimUnroll, dims),
                                                             "unrollDim");
      sourceStream << this->indent[1] << "}" << std::endl;

      if (dims % maxDimUnroll != 0) {
        sourceStream << this->unrolledBasisFunctionEvalulation(
            dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "");
      }
    } else {
      sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");
    }

    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[1] << "myResult_" << i << " += curSupport_" << i << ";"
                   << std::endl;
    }
    sourceStream << this->indent[0] << "}" << std::endl;
    sourceStream << std::endl;

    if (useLocalMemory) {
      sourceStream << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << this->indent[0] << "}" << std::endl;
    }

    sourceStream << std::endl;
    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[0] << "ptrResult[(globalSize * " << i
                   << ") + globalIdx] = myResult_" << i << ";" << std::endl;
    }
    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("StreamingOCLMultiPlatform_mult.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};
}  // namespace StreamingOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
