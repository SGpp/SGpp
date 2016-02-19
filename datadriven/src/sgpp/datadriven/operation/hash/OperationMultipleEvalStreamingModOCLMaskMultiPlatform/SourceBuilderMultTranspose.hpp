// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <fstream>
#include <memory>
#include <string>

#include "sgpp/base/exception/operation_exception.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/base/opencl/OCLDevice.hpp"
#include "sgpp/base/opencl/KernelSourceBuilderBase.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingModOCLMaskMultiPlatform {

template <typename T>
class SourceBuilderMultTranspose : public base::KernelSourceBuilderBase<T> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;

 public:
  SourceBuilderMultTranspose(std::shared_ptr<base::OCLDevice> device,
                             json::Node &kernelConfiguration, size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
  }

  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("streamingModOCLMaskMP_multTranspose.cl");
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
    sourceStream << "void multTransOCLMask(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrMask," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrOffset,"
                 << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrSource,"
                 << std::endl;
    sourceStream << "           __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    sourceStream << "           uint start_data," << std::endl;
    sourceStream << "           uint end_data)" << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << "   int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << "   int localIdx = get_local_id(0);" << std::endl;
    sourceStream << "   uint rangeData = end_data - start_data;" << std::endl;

    sourceStream << std::endl;
    sourceStream << "   " << this->floatType()
                 << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl
                 << std::endl;
    sourceStream << "   " << this->floatType() << " myResult = ptrResult[globalIdx];" << std::endl
                 << std::endl;
    if (useLocalMemory) {
      sourceStream << "   __local " << this->floatType() << " locData[" << dims * localWorkgroupSize
                   << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locSource[" << localWorkgroupSize
                   << "];" << std::endl
                   << std::endl;
    }

    for (size_t d = 0; d < dims; d++) {
      sourceStream << " " << this->floatType() << " level_" << d << " = ptrLevel[(globalIdx*"
                   << dims << ")+" << d << "];" << std::endl;
      sourceStream << " " << this->floatType() << " index_" << d << " = ptrIndex[(globalIdx*"
                   << dims << ")+" << d << "];" << std::endl;
      sourceStream << " " << this->floatType() << " mask_" << d << " = ptrMask[(globalIdx*" << dims
                   << ")+" << d << "];" << std::endl;
      sourceStream << " " << this->floatType() << " offset_" << d << " = ptrOffset[(globalIdx*"
                   << dims << ")+" << d << "];" << std::endl;
    }

    sourceStream << std::endl;
    sourceStream << "   // Iterate over all grid points" << std::endl;
    if (useLocalMemory) {
      sourceStream << " for(int i = start_data; i < end_data; i+=" << localWorkgroupSize << ")"
                   << std::endl;
      sourceStream << "   {" << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << "     locData[(" << d << "*" << localWorkgroupSize
                     << ")+(localIdx)] = ptrData[(" << d << "*rangeData)+(localIdx+i)];"
                     << std::endl;
      }

      sourceStream << "       locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
      sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << std::endl;
      sourceStream << "       for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
      sourceStream << "       {" << std::endl;

      sourceStream << "           curSupport = locSource[k];" << std::endl
                   << std::endl;
    } else {
      sourceStream << "   for(int k = start_data; k < end_data; k++)" << std::endl;
      sourceStream << "   {" << std::endl;
      sourceStream << "     curSupport = ptrSource[k];" << std::endl
                   << std::endl;
    }

    for (size_t d = 0; d < dims; d++) {
      if (useLocalMemory) {
        sourceStream << "         eval = ((level_" << d << ") * (locData[(" << d << "*"
                     << localWorkgroupSize << ")+k]));" << std::endl;
      } else {
        sourceStream << "         eval = ((level_" << d << ") * (ptrData[(" << d
                     << "*rangeData)+k]));" << std::endl;
      }
      sourceStream << "         index_calc = eval - (index_" << d << ");" << std::endl;
      sourceStream << "         abs = as_" << this->floatType() << "(as_" << this->intType()
                   << "(index_calc) | as_" << this->intType() << "(mask_" << d << "));"
                   << std::endl;
      sourceStream << "         last = offset_" << d << " + abs;" << std::endl;
      sourceStream << "         localSupport = fmax(last, 0.0" << this->constSuffix() << ");"
                   << std::endl;
      sourceStream << "         curSupport *= localSupport;" << std::endl;
    }

    sourceStream << std::endl
                 << "      myResult += curSupport;" << std::endl;
    sourceStream << "       }" << std::endl
                 << std::endl;

    if (useLocalMemory) {
      sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << "   }" << std::endl;
    }

    sourceStream << "   ptrResult[globalIdx] = myResult;" << std::endl;
    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("streamingModOCLMaskMP_multTranspose.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};

}  // namespace StreamingModOCLMaskMultiPlatform
}  // namespace datadriven
}  // namespace SGPP
