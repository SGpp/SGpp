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

/// OpenCL source builder for density based graph pruning
template<typename real_type>
class SourceBuilderPruneGraph: public base::KernelSourceBuilderBase<real_type> {
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

 public:
  SourceBuilderPruneGraph(json::Node &kernelConfiguration, size_t dims) :
      kernelConfiguration(kernelConfiguration), dims(dims) {
  }

  /// Generates the whole opencl kernel code used for the pruning of a graph
  std::string generateSource(size_t dimensions, size_t gridSize, size_t k, real_type treshold) {
    if (kernelConfiguration.contains("REUSE_SOURCE")) {
      if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
        return this->reuseSource("DensityOCLMultiPlatform_prune_graph.cl");
      }
    }

    std::stringstream sourceStream;

    if (this->floatType().compare("double") == 0) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable"
                   << std::endl << std::endl;
    }

    sourceStream << "" << std::endl
                 << "" << this->floatType() << " get_u(__private const "
                 << this->floatType() << " grenze,__private const int index," << std::endl
                 << "__private const int level)"
                 << " {" << std::endl
                 << this->indent[0] << "private " << this->floatType()
                 << " ret = (1 << level);" << std::endl
                 << this->indent[0] << "ret*=grenze;" << std::endl
                 << this->indent[0] << "ret-=index;" << std::endl
                 << this->indent[0] << "ret = fabs(ret);" << std::endl
                 << this->indent[0] << "ret=1-ret;" << std::endl
                 << this->indent[0] << "if (ret<0.0)" << std::endl
                 << this->indent[1] << "ret=0.0;" << std::endl
                 << this->indent[0] << "return ret;" << std::endl
                 << "}" << std::endl
                 << "" << std::endl
                 << "void kernel removeEdges(__global int *nodes,"
                 << "__global const int *starting_points,__global const " << this->floatType()
                 << " *data," << std::endl
                 << this->indent[0] << "__global const " << this->floatType()
                 << " *alphas, int startid) {" << std::endl
                 << this->indent[0] << "size_t index = get_global_id(0);" << std::endl
                 << this->indent[0] << "size_t global_index = startid + get_global_id(0);"
                 << std::endl
                 << this->indent[0] << "__private " << this->floatType()
                 << " endwert=0.0;" << std::endl
                 << this->indent[0] << "__private " << this->floatType()
                 << " wert=1.0;" << std::endl
                 << this->indent[0] << "for (int i = 0; i <  " << k << " ; i++)" << std::endl
                 << this->indent[0] << "{" << std::endl
                 << this->indent[1] << "//Calculate density" << std::endl
                 << this->indent[1] << "endwert=0;" << std::endl
                 << this->indent[1] << "int nachbar=nodes[index* " << k << " +i];" << std::endl
                 << this->indent[1] << "for (int gridpoint=0;gridpoint< " << gridSize
                 << " ;gridpoint++)" << std::endl
                 << this->indent[1] << "{" << std::endl
                 << this->indent[2] << "wert=1;" << std::endl
                 << this->indent[2] << "for (int dimension=0;dimension< " << dimensions
                 << ";dimension++)" << std::endl
                 << this->indent[2] << "{" << std::endl
                 << this->indent[3] << "" << this->floatType() << " dimension_point=0;"
                 << std::endl
                 << this->indent[4] << "dimension_point=data[dimension+nachbar* "
                 << dimensions << "]+" << std::endl
                 << this->indent[4] << "(data[global_index* " << dimensions
                 << "+dimension]-data[dimension+nachbar* " << dimensions << "])*0.5;"
                 << std::endl
                 << this->indent[3] << "wert*=get_u(dimension_point,"
                 << "starting_points[gridpoint*2* " << dimensions << "+2*dimension],"
                 << std::endl
                 << this->indent[3] << "starting_points[gridpoint*2* " << dimensions
                 << "+2*dimension+1]);" << std::endl
                 << this->indent[2] << "}" << std::endl
                 << this->indent[2] << "endwert+=wert*alphas[gridpoint];" << std::endl
                 << this->indent[1] << "}" << std::endl
                 << this->indent[1] << "if (endwert< " << treshold << " )" << std::endl
                 << this->indent[1] << "{" << std::endl
                 << this->indent[2] << "nodes[ " << k << " *index + i] = -2;" << std::endl
                 << this->indent[1] << "}" << std::endl
                 << this->indent[0] << "}" << std::endl
                 << this->indent[0] << "endwert=0;" << std::endl
                 << this->indent[0] << "for (int gridpoint=0;gridpoint< " << gridSize
                 << " ;gridpoint++)" << std::endl
                 << this->indent[0] << "{" << std::endl
                 << this->indent[1] << "wert=1;" << std::endl
                 << this->indent[1] << "for (int dimension=0;dimension< " << dimensions
                 << ";dimension++)" << std::endl
                 << this->indent[1] << "{" << std::endl
                 << this->indent[2] << "wert*=get_u(data[global_index* " << dimensions
                 << "+dimension],starting_points[gridpoint*2* " << dimensions
                 << "+2*dimension]," << std::endl
                 << this->indent[3] << "starting_points[gridpoint*2* " << dimensions
                 << "+2*dimension+1]);" << std::endl
                 << this->indent[1] << "}" << std::endl
                 << this->indent[1] << "endwert+=wert*alphas[gridpoint];" << std::endl
                 << this->indent[0] << "}" << std::endl
                 << this->indent[0] << "if (endwert< " << treshold << " )" << std::endl
                 << this->indent[0] << "{" << std::endl
                 << this->indent[1] << "for (int i = 0; i <  " << k << " ; i++)" << std::endl
                 << this->indent[1] << "{" << std::endl
                 << this->indent[2] << "nodes[ " << k << " *index + i] = -1;" << std::endl
                 << this->indent[1] << "}" << std::endl
                 << this->indent[0] << "}" << std::endl
                 << "}" << std::endl;
    if (kernelConfiguration.contains("WRITE_SOURCE")) {
      if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
        this->writeSource("DensityOCLMultiPlatform_prune_graph.cl", sourceStream.str());
      }
    }
    return sourceStream.str();
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
