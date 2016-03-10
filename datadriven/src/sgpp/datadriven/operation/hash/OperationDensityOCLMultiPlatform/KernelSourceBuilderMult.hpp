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

template<typename real_type>
class SourceBuilderMult: public base::KernelSourceBuilderBase<real_type> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
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
  SourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration,
                    size_t dims) :
      device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
  }

  std::string generateSource(size_t dimensions, size_t gridsize) {
    std::stringstream sourceStream;

    if (this->floatType().compare("double") == 0) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl
                   << std::endl;
    }

    sourceStream << "void kernel multdensity(global const int *starting_points,global const "
                 << this->floatType() << " *alpha, global " << this->floatType()
                 << " *result,const " << this->floatType() << " lambda, int startid)"
                 << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << this->indent[0] << "int gridindex = startid + get_global_id(0);"
                 << std::endl;
    sourceStream << this->indent[0] << this->floatType() << " gesamtint = 0.0;" << std::endl;
    sourceStream << this->indent[0] << "for(private int i = 0;i< " << gridsize
                 << ";i++) {" << std::endl;
    sourceStream << this->indent[1] << this->floatType() << " zellenintegral = 1.0;"
                 << std::endl;
    sourceStream << this->indent[1] << "for(private int dim = 0;dim< " << dimensions
                 << ";dim++) {" << std::endl;
    sourceStream << this->indent[2] << "int index = starting_points[gridindex* "
                 << dimensions << "*2+2*dim];" << std::endl;
    sourceStream << this->indent[2] << "int level = starting_points[gridindex* "
                 << dimensions << "*2+2*dim+1];" << std::endl;
    sourceStream << this->indent[2] << "int index2 = starting_points[i* " << dimensions
                 << "*2+2*dim];" << std::endl;
    sourceStream << this->indent[2] << "int level2 = starting_points[i* " << dimensions
                 << "*2+2*dim+1];" << std::endl;
    sourceStream << this->indent[2] << "if(starting_points[gridindex* " << dimensions
                 << "*2+2*dim+1]>starting_points[i* " << dimensions
                 << "*2+2*dim+1]) {" << std::endl;
    sourceStream << this->indent[3] << "index = starting_points[i* "
                 << dimensions << "*2+2*dim];" << std::endl;
    sourceStream << this->indent[3] << "level = starting_points[i* "
                 << dimensions << "*2+2*dim+1];" << std::endl;
    sourceStream << this->indent[3] << "index2 = starting_points[gridindex* "
                 << dimensions << "*2+2*dim];" << std::endl;
    sourceStream << this->indent[3] << "level2 = starting_points[gridindex* "
                 << dimensions << "*2+2*dim+1];" << std::endl;
    sourceStream << this->indent[2] << "}" << std::endl;
    sourceStream << this->indent[2] << "int teiler = (1 << level2);" << std::endl;
    sourceStream << this->indent[2] << this->floatType() << " h = 1.0 / teiler;" << std::endl;
    sourceStream << this->indent[2] << this->floatType() << " grenze1 = h*(index2-1);"
                 << std::endl;
    sourceStream << this->indent[2] << this->floatType() << " grenze2 = h*(index2+1);"
                 << std::endl;
    sourceStream << this->indent[2] << "int u= (1 << level);" << std::endl;
    sourceStream << this->indent[2] << this->floatType() << " uright = u*grenze2-index;"
                 << std::endl;
    sourceStream << this->indent[2] << this->floatType() << " uleft = u*grenze1-index;"
                 << std::endl;
    sourceStream << this->indent[2] << "if(uleft<0)" << std::endl;
    sourceStream << this->indent[3] << "uleft *= -1;" << std::endl;
    sourceStream << this->indent[2] << "uleft = 1-uleft;" << std::endl;
    sourceStream << this->indent[2] << "if(uleft<0)" << std::endl;
    sourceStream << this->indent[3] << "uleft = 0;" << std::endl;
    sourceStream << this->indent[2] << "if(uright<0)" << std::endl;
    sourceStream << this->indent[3] << "uright *= -1;" << std::endl;
    sourceStream << this->indent[2] << "uright = 1-uright;" << std::endl;
    sourceStream << this->indent[2] << "if(uright<0)" << std::endl;
    sourceStream << this->indent[3] << "uright = 0;" << std::endl;
    sourceStream << this->indent[2] << "if(starting_points[i* " << dimensions
                 << "*2+2*dim+1] == starting_points[gridindex* "
                 << dimensions << "*2+2*dim+1]) {" << std::endl;
    sourceStream <<  this->indent[3] <<"zellenintegral *= 2.0/3.0*h;" << std::endl;
    sourceStream << this->indent[3] << "if(starting_points[i* " << dimensions
                 << "*2+2*dim] != starting_points[gridindex* "
                 << dimensions << "*2+2*dim])" << std::endl;
    sourceStream << this->indent[4] << "zellenintegral = 0.0;" << std::endl;
    sourceStream << this->indent[2] << "}" << std::endl;
    sourceStream << this->indent[2] << "else" << std::endl;
    sourceStream << this->indent[3] << "zellenintegral *= h/2.0*(uleft+uright);" << std::endl;
    sourceStream << this->indent[1] << "}" << std::endl;
    sourceStream << this->indent[1] << "gesamtint += zellenintegral*alpha[i];" << std::endl;
    sourceStream << this->indent[0] << "}" << std::endl;
    sourceStream << this->indent[0] << "result[get_global_id(0)] = gesamtint;" << std::endl;
    sourceStream << this->indent[0] << "result[get_global_id(0)] += alpha[gridindex]*"
                 << "lambda;" << std::endl;
    sourceStream << "}" << std::endl;

    if (kernelConfiguration.contains("WRITE_SOURCE")) {
      if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
        this->writeSource("DensityOCLMultiPlatform_mult.cl", sourceStream.str());
      }
    }
    return sourceStream.str();
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
