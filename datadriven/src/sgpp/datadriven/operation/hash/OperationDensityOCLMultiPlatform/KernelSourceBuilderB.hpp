// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

#include <fstream>
#include <sstream>
#include <string>
namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// OpenCL source builder for density right hand side vector
template<typename real_type>
class SourceBuilderB: public base::KernelSourceBuilderBase<real_type> {
 private:
  /// OpenCL configuration containing the building flags
  json::Node &kernelConfiguration;
  /// Dimensions of grid
  size_t dims;
  /// Used workgroupsize for opencl kernel execution
  size_t localWorkgroupSize;
  /// Using local memory?
  bool useLocalMemory;
  size_t dataBlockSize;
  size_t transGridBlockSize;
  uint64_t maxDimUnroll;

  std::string getData(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      output << "data_" << dataBlockingIndex << "[" << dim << "]";
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      output << "data_" << dataBlockingIndex << "_" << dim;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim
             << ") + " << dataBlockingIndex << "]";
    } else {
      std::string error("OCL Error: Illegal value for parameter \"KERNEL_STORE_DATA\"");
      throw new base::operation_exception(error.c_str());
    }
    return output.str();
  }

 public:
  SourceBuilderB(json::Node &kernelConfiguration,
                 size_t dims) :
      kernelConfiguration(kernelConfiguration), dims(dims) {
  }

  /// Generates the opencl source code for the density right hand side vector
  std::string generateSource(size_t dimensions, size_t datapoints) {
    if (kernelConfiguration.contains("REUSE_SOURCE")) {
      if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
        return this->reuseSource("DensityOCLMultiPlatform_mult.cl");
      }
    }

    std::stringstream sourceStream;

    if (this->floatType().compare("double") == 0) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl
                   << std::endl;
    }

    sourceStream << "void kernel cscheme(global const int* starting_points," << std::endl
                 <<"global const " << this->floatType() << "* data_points,global "
                 << this->floatType() << "* C, private int startid) {" << std::endl
                 << this->indent[0] << "C[get_global_id(0)]=0.0;" << std::endl
                 << this->indent[0] << "private " << this->floatType() << " value=1;"
                 << std::endl
                 << this->indent[0] << "private " << this->floatType() << " wert=1.0;"
                 << std::endl
                 << this->indent[0] << "for(unsigned int ds=0;ds< " << datapoints << ";ds++)"
                 << std::endl
                 << this->indent[0] << "{" << std::endl
                 << this->indent[1] << "value=1;"
                 << std::endl
                 << this->indent[1] << "for(private int d=0;d< " << dimensions << ";d++)"
                 << std::endl
                 << this->indent[1] << "{" << std::endl
                 << this->indent[2] <<"wert = (1 << starting_points[(startid + "
                 << " get_global_id(0))*2* " << dimensions << "+2*d+1]);" << std::endl
                 << this->indent[2] << "wert*=data_points[ds* " << dimensions << "+d];"
                 << std::endl
                 << this->indent[2] << "wert-=starting_points[(startid + get_global_id(0))*2* "
                 << dimensions << "+2*d];" << std::endl
                 << this->indent[3] << "wert=fabs(wert);" << std::endl
                 << this->indent[2] << "wert=1-wert;" << std::endl
                 << this->indent[2] << "if(wert<0)" << std::endl
                 << this->indent[3] << "wert=0;" << std::endl
                 << this->indent[2] << "value*=wert;" << std::endl
                 << this->indent[1] << "}" << std::endl
                 << this->indent[1] << "C[get_global_id(0)]+=value;"
                 << std::endl
                 << this->indent[0] << "}" << std::endl
                 << this->indent[1] << "C[get_global_id(0)]/=" << datapoints << ";"
                 << std::endl
                 <<"}" << std::endl;
    if (kernelConfiguration.contains("WRITE_SOURCE")) {
      if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
        this->writeSource("DensityOCLMultiPlatform_rhs.cl", sourceStream.str());
      }
    }
    return sourceStream.str();
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
