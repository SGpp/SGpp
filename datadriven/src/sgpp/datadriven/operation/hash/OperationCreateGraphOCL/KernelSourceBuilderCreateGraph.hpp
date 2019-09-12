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

/// OpenCL source builder for k nearest neighbors graph creation
template<typename real_type>
class SourceBuilderCreateGraph: public base::KernelSourceBuilderBase<real_type> {
 private:
  /// OpenCL configuration containing the building flags
  json::Node &kernelConfiguration;
  /// Dimensions of the used dataset
  size_t dims;
  /// Used workgroupsize for opencl kernel execution
  size_t localWorkgroupSize;
  /// Using local memory?
  bool useLocalMemory;
  size_t dataBlockSize;
  size_t transGridBlockSize;
  uint64_t maxDimUnroll;
  /// Use select statements instead of if branches? Configuration parameter is USE_SELECT
  bool use_select;
  /// Create a approximation instead of an exact graph? Configuration parameter is USE_APPROX
  bool use_approx;
  /// Number of bins used for the approximation if use_approx is true
  size_t approxRegCount;

  /// Writes the source code for initialization of the k neighbor registers
  std::string init_k_registers(size_t k) {
    std::stringstream output;
    output << this->indent[0] << "__private int k_reg["<< k << "];" << std::endl;
    output << this->indent[0] << "__private " << this->floatType() << " k_dists["
           << k << "];" << std::endl;
    output << this->indent[0] << "for (int i = 0; i < " << k << "; i++)" << std::endl;
    output << this->indent[1] << "k_dists[i] = 4.0;" << std::endl;
    /*for (size_t i = 0; i < k; i++) {
      output << this->indent[0] << "int k_register" << i << " = " << i << "; " << std::endl;
      }
      for (size_t i = 0; i < k; i++) {
      output << this->indent[0] << this->floatType()
      << " dist_k" << i << " = 4.0;" << std::endl;
      }*/
    return output.str();
  }
  /// Writes the source code for finding the current maximum of all neighbors
  std::string find_max_index(size_t k, bool unroll) {
    std::stringstream output;
    if (!unroll) {
      output << this->indent[2] << "for (int j = 1; j < "
             << k << "; j++) {" << std::endl
             << this->indent[3] << "if (k_dists[maxindex] < k_dists[j])" << std::endl
             << this->indent[4] << "maxindex = j;" << std::endl
             << this->indent[2] << "}" << std::endl;
    } else {
      for (size_t i = 1; i < k; i++) {
        if (use_select) {
          output << this->indent[2] << "maxindex = select(maxindex, " << i
                 << ", k_dists[maxindex] < k_dists["
                 << i << "]);" << std::endl;
        } else {
          output << this->indent[2] << "if (k_dists[maxindex] < k_dists[" << i << "])" << std::endl;
          output << this->indent[3] << "maxindex  = " << i << "; " << std::endl;
        }
      }
    }
    // Enables vectorization but slows kernel down (longer ifs...)
    /*for (size_t i = 0; i < k; i++) {
      output << this->indent[2] << "if (maxdist < k_dists[" << i << "]) {" << std::endl;
      output << this->indent[3] << "maxindex  = " << i << "; " << std::endl;
      output << this->indent[3] << "maxdist  = k_dists["  << i << "]; }" << std::endl;
      }*/
    return output.str();
  }
  /// Writes the source code for copying the current datapoint to private memory
  std::string save_from_global_to_private(size_t dimensions) {
    std::stringstream output;
    output << this->indent[0] << "__private " << this->floatType() << " datapoint["
           << dimensions << "];" << std::endl;
    for (size_t i = 0; i < dimensions; i++) {
      output << this->indent[1] << "datapoint[" << i << "] = data[global_index * "
             << dimensions << " + " << i << "];" << std::endl;
    }
    return output.str();
  }
  /// Writes the source code for copying the results back into the global memory
  std::string copy_k_registers_into_global(size_t k) {
    std::stringstream output;
    for (size_t i = 0; i < k; i++) {
      output << this->indent[0] << "neighbors[chunk_index * "<< k
             << " + " << i << "] = k_reg[" << i << "];" << std::endl;
    }
    return output.str();
  }

 public:
  SourceBuilderCreateGraph(json::Node &kernelConfiguration, size_t dims) :
      kernelConfiguration(kernelConfiguration), dims(dims), use_select(false),
      use_approx(false) {
    localWorkgroupSize = 128;
    if (kernelConfiguration.contains("LOCAL_SIZE"))
      localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    approxRegCount = localWorkgroupSize;
    if (kernelConfiguration.contains("KERNEL_USE_LOCAL_MEMORY"))
      useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    if (kernelConfiguration.contains("USE_SELECT")) {
      if (kernelConfiguration["USE_SELECT"].getBool()) {
        use_select = true;
      }
    }
    if (kernelConfiguration.contains("USE_APPROX")) {
      if (kernelConfiguration["USE_APPROX"].getBool()) {
        use_approx = true;
      }
    }
    if (kernelConfiguration.contains("APPROX_REG_COUNT"))
      approxRegCount = kernelConfiguration["APPROX_REG_COUNT"].getUInt();
  }

  /// Generates the whole opencl kernel code for the creation of a k nearest neighbor graph
  std::string generateSource(size_t dimensions, size_t k, size_t dataSize, size_t problemsize) {
    if (kernelConfiguration.contains("REUSE_SOURCE")) {
      if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
        return this->reuseSource("DensityOCLMultiPlatform_create_graph.cl");
      }
    }

    std::stringstream sourceStream;

    if (this->floatType().compare("double") == 0) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable"
                   << std::endl << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl
                 << "void connectNeighbors(__global const " << this->floatType()
                 << " *data, __global int *neighbors, const int startid)"
                 << std::endl
                 << "{" << std::endl
                 << this->indent[0] << "__private int global_index = startid + get_global_id(0);"
                 << this->indent[0] << "__private int local_id = get_local_id(0);"
                 << std::endl
                 << this->indent[0] << "__private int chunk_index = get_global_id(0);" << std::endl
                 << this->indent[0] << "__private int maxindex = 0;" << std::endl;
    if (!use_approx)
      sourceStream << init_k_registers(k) << std::endl;
    sourceStream << save_from_global_to_private(dimensions)
                 << this->indent[0] << "__private " << this->floatType() << " dist = 0.0;"
                 << std::endl;
    if (use_approx) {
      sourceStream << this->indent[0] << "__local " << this->floatType() << " data_local["
                   << localWorkgroupSize * dimensions  << "];" << std::endl
                   << this->indent[0] << "__private " << this->floatType() << " dist_reg["
                   << approxRegCount  << "];" << std::endl
                   << this->indent[0] << "__private int  index_reg["
                   << approxRegCount  << "];" << std::endl
                   << this->indent[1] << "for (int i = 0 ; i < " << approxRegCount
                   << "; i++) " << std::endl
                   << this->indent[1] << "dist_reg[i] = 4.0;" << std::endl
                   <<  this->indent[0] << "for (int group = 0; group < "
                   << problemsize / localWorkgroupSize << "; group++) {" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int j = 0; j <     " << dimensions
                   << " ; j++) " << std::endl
                   << this->indent[2] << "data_local[local_id * " << dimensions
                   << " + j] = data[group * " << localWorkgroupSize * dimensions
                   << "  + local_id * " <<  dimensions << " + j];" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      if (localWorkgroupSize == approxRegCount) {
        sourceStream << this->indent[1] << "for (int i = 0 ; i < " << localWorkgroupSize
                     << "; i++) {" << std::endl;
      } else {
        sourceStream << this->indent[2] << "int chunkindex = 0;" << std::endl;
        sourceStream << this->indent[1] << "for (int chunk = 0 ; chunk < "
                     << localWorkgroupSize / approxRegCount << "; chunk++) {" << std::endl;
        sourceStream << this->indent[1] << "for (int i = 0 ; i < " << approxRegCount
                     << "; i++, chunkindex++) {" << std::endl;
      }
      sourceStream << this->indent[2] << "dist = 0.0;" << std::endl
                   << this->indent[2] << "for (int j = 0; j <     " << dimensions
                   << " ; j++) {" << std::endl;
      if (localWorkgroupSize != approxRegCount) {
        sourceStream << this->indent[3] << "dist += (datapoint[j] - data_local[j + (chunkindex) * "
                     << dimensions << " ])" << std::endl
                     << this->indent[3] << "* (datapoint[j] - data_local[j + (chunkindex)* "
                     << dimensions << " ]);" << std::endl;
        sourceStream << this->indent[2] << "}" << std::endl
                     << this->indent[2] << "if (dist < dist_reg[i] && chunkindex + group * "
                     << localWorkgroupSize << " != global_index) {" << std::endl
                     << this->indent[3] << "dist_reg[i] = dist;" << std::endl
                     << this->indent[3] << "index_reg[i] = group * " << localWorkgroupSize
                     << " + chunkindex;" << std::endl;
      } else {
        sourceStream << this->indent[3] << "dist += (datapoint[j] - data_local[j + i * "
                     << dimensions << " ])" << std::endl
                     << this->indent[3] << "* (datapoint[j] - data_local[j + i* "
                     << dimensions << " ]);" << std::endl;
        sourceStream << this->indent[2] << "}" << std::endl
                     << this->indent[2] << "if (dist < dist_reg[i] && i + group * "
                     << localWorkgroupSize << " != global_index) {" << std::endl
                     << this->indent[3] << "dist_reg[i] = dist;" << std::endl
                     << this->indent[3] << "index_reg[i] = group;" << std::endl;
      }
      sourceStream << this->indent[2] << "}" << std::endl
                   << this->indent[1] << "}" << std::endl;
      if (localWorkgroupSize != approxRegCount)
        sourceStream << this->indent[1] << "}" << std::endl;
      sourceStream << this->indent[0] << "}" << std::endl
                   << this->indent[0] << "for (int neighbor = 0 ; neighbor < " << k
                   << "; neighbor++) {" << std::endl
                   << this->indent[1] << "maxindex = 0;" << std::endl
                   << this->indent[1] << "for (int i = 1 ; i < " << approxRegCount
                   << "; i++) {" << std::endl
                   << this->indent[2] << "if (dist_reg[i] < dist_reg[maxindex])" << std::endl
                   << this->indent[3] << "maxindex = i;" << std::endl
                   << this->indent[1] << "}" << std::endl;

      if (localWorkgroupSize == approxRegCount) {
        sourceStream << this->indent[1] << "neighbors[chunk_index * "<< k
                     << " + neighbor] = maxindex + index_reg[maxindex] * "
                     << localWorkgroupSize << ";" << std::endl;
      } else {
        sourceStream << this->indent[1] << "neighbors[chunk_index * "<< k
                     << " + neighbor] = index_reg[maxindex];" << std::endl;
      }
      sourceStream << this->indent[1] << "dist_reg[maxindex] = 4.0;" << std::endl;
      sourceStream << this->indent[0] << "}" << std::endl;
    } else if (useLocalMemory) {
      sourceStream << this->indent[0] << "__local " << this->floatType() << " data_local["
                   << localWorkgroupSize * dimensions << "];" << std::endl
                   <<  this->indent[0] << "for (int group = 0; group < "
                   << problemsize / localWorkgroupSize << "; group++) {" << std::endl
                   << this->indent[1] << "for (int j = 0; j <     " << dimensions
                   << " ; j++) " << std::endl
                   << this->indent[2] << "data_local[local_id * " << dimensions
                   << " + j] = data[group * " << localWorkgroupSize * dimensions
                   << "  + local_id * " <<  dimensions << " + j];" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int i = 0 ; i < " << localWorkgroupSize
                   << "; i++) {" << std::endl
                   << this->indent[2] << "dist = 0.0;" << std::endl
                   << this->indent[2] << "for (int j = 0; j <     " << dimensions
                   << " ; j++) {" << std::endl
                   << this->indent[3] << "dist += (datapoint[j] - data_local[j + i * "
                   << dimensions << " ])" << std::endl
                   << this->indent[3] << "* (datapoint[j] - data_local[j + i* "
                   << dimensions << " ]);" << std::endl
                   << this->indent[2] << "}" << std::endl
                   << this->indent[2] << "maxindex = 0;" << std::endl
                   << find_max_index(k, true)
                   << this->indent[2] << "if (dist < k_dists[maxindex] && i + group * "
                   << localWorkgroupSize << " != global_index) {" << std::endl
                   << this->indent[3] << "k_reg[maxindex] = i + group * " << localWorkgroupSize
                   << ";" << std::endl
                   << this->indent[3] << "k_dists[maxindex] = dist;" << std::endl
                   << this->indent[2] << "}" << std::endl
                   << this->indent[1] << "}" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[0] << "}" << std::endl
                   << copy_k_registers_into_global(k);
    } else {
      sourceStream << this->indent[0] << "for (int i = 0; i <    " << dataSize
                   << "; i++) {" << std::endl
                   << this->indent[1] << "if (i != global_index) {" << std::endl
                   << "//get distance to current point" << std::endl
                   << this->indent[2] << "dist = 0.0;" << std::endl
                   << this->indent[2] << "for (int j = 0; j <     " << dimensions
                   << " ; j++) {" << std::endl
                   << this->indent[3] << "dist += (datapoint[j] - data[j + i* "
                   << dimensions << " ])" << std::endl
                   << this->indent[3] << "* (datapoint[j] - data[j + i* "
                   << dimensions << " ]);" << std::endl
                   << this->indent[2] << "}" << std::endl
                   << this->indent[2] << "int maxindex = 0;" << std::endl
                   << find_max_index(k, true)
                   << this->indent[2] << "if (dist < k_dists[maxindex]) {" << std::endl
                   << this->indent[3] << "k_reg[maxindex] = i;" << std::endl
                   << this->indent[3] << "k_dists[maxindex] = dist;" << std::endl
                   << this->indent[2] << "}" << std::endl
                   << this->indent[1] << "}" << std::endl
                   << this->indent[0] << "}" << std::endl
                   << copy_k_registers_into_global(k);
    }
    sourceStream << "}" << std::endl;
    if (kernelConfiguration.contains("WRITE_SOURCE")) {
      if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
        this->writeSource("DensityOCLMultiPlatform_create_graph.cl", sourceStream.str());
      }
    }
    return sourceStream.str();
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
