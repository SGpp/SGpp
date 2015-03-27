/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>

//#include "StreamingOCLParameters.hpp"

namespace SGPP {
namespace datadriven {

template<typename T> struct getType {
};

template<> struct getType<float> {
  static std::string asString() {
    return "float";
  }
  static std::string constSuffix() {
    return "f";
  }
  static std::string intAsString() {
    return "uint";
  }
};

template<> struct getType<double> {
  static std::string asString() {
    return "double";
  }
  static std::string constSuffix() {
    return "";
  }
  static std::string intAsString() {
    return "ulong";
  }
};

template<typename real_type>
class OCLKernelSourceBuilder {
private:
  base::ConfigurationParameters parameters;
  real_type blubb;
public:
  OCLKernelSourceBuilder(base::ConfigurationParameters parameters) :
      parameters(parameters) {

  }

  std::string generateSourceMult(size_t dims) {

    size_t localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    bool useLocalMemory = this->parameters.getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY");
    uint64_t maxDimUnroll = this->parameters.getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL");

    std::stringstream stream_program_src;

    if (getType<real_type>::asString() == "double") {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    stream_program_src << "__kernel" << std::endl;
    stream_program_src << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    stream_program_src << "void multOCL(__global const " << getType<real_type>::asString() << "* ptrLevel,"
        << std::endl;
    stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrIndex," << std::endl;
    stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrData," << std::endl;
    stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrAlpha," << std::endl;
    stream_program_src << "           __global       " << getType<real_type>::asString() << "* ptrResult," << std::endl;
    stream_program_src << "           uint resultSize," << std::endl;
    stream_program_src << "           uint start_grid," << std::endl;
    stream_program_src << "           uint end_grid) " << std::endl;
    stream_program_src << "{" << std::endl;
    stream_program_src << " int globalIdx = get_global_id(0);" << std::endl;
    stream_program_src << " int localIdx = get_local_id(0);" << std::endl;

    stream_program_src << std::endl;

    if (useLocalMemory) {
      stream_program_src << " __local " << getType<real_type>::asString() << " locLevel[" << dims * localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << " __local " << getType<real_type>::asString() << " locIndex[" << dims * localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << " __local " << getType<real_type>::asString() << " locAlpha[" << localWorkgroupSize << "];"
          << std::endl;
      stream_program_src << std::endl;
    }
    stream_program_src << " " << getType<real_type>::asString()
        << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
    stream_program_src << " " << getType<real_type>::asString() << " myResult = 0.0;" << std::endl << std::endl;
    stream_program_src << " // Create registers for the data" << std::endl;

    if (dims > maxDimUnroll) {
      stream_program_src << " " << getType<real_type>::asString() << " data[" << dims << "];" << std::endl;
      stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
      stream_program_src << "  data[d] = ptrData[globalIdx+(resultSize * d)];" << std::endl;
      stream_program_src << "}" << std::endl;
    } else {
      for (size_t d = 0; d < dims; d++) {
        stream_program_src << " " << getType<real_type>::asString() << " data_" << d
            << " = ptrData[globalIdx+(resultSize*" << d << ")];" << std::endl;
      }
    }

    stream_program_src << std::endl;

    if (useLocalMemory) {
      stream_program_src << " // Iterate over all grid points (fast ones, with cache)" << std::endl;
      stream_program_src << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
      stream_program_src << " uint fastChunkSizeGrid = (chunkSizeGrid / " << localWorkgroupSize << ") * "
          << localWorkgroupSize << ";" << std::endl;
      stream_program_src << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+=" << localWorkgroupSize
          << ")" << std::endl;
      stream_program_src << " {" << std::endl;

      if (dims > maxDimUnroll) {
        stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "   locLevel[(localIdx*" << dims << ")+ d] = ptrLevel[((j+localIdx)*" << dims << ")+ d];"
            << std::endl;
        stream_program_src << "   locIndex[(localIdx*" << dims << ")+ d] = ptrIndex[((j+localIdx)*" << dims << ")+ d];"
            << std::endl;
        stream_program_src << "}" << std::endl;
      } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "   locLevel[(localIdx*" << dims << ")+" << d << "] = ptrLevel[((j+localIdx)*" << dims
              << ")+" << d << "];" << std::endl;
          stream_program_src << "   locIndex[(localIdx*" << dims << ")+" << d << "] = ptrIndex[((j+localIdx)*" << dims
              << ")+" << d << "];" << std::endl;
        }
      }

      stream_program_src << "   locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
      stream_program_src << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream_program_src << std::endl;
      stream_program_src << "   for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
      stream_program_src << "   {" << std::endl;
      stream_program_src << "     curSupport = locAlpha[k];" << std::endl << std::endl;

      if (dims > maxDimUnroll) {
        stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "     eval = ((locLevel[(k*" << dims << ")+ d ]) * (data[d]));" << std::endl;
        stream_program_src << "     index_calc = eval - (locIndex[(k*" << dims << ")+ d]);" << std::endl;
        stream_program_src << "     abs = fabs(index_calc);" << std::endl;
        stream_program_src << "     last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
        stream_program_src << "     localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
            << std::endl;
        stream_program_src << "     curSupport *= localSupport;" << std::endl << std::endl;
        stream_program_src << "}" << std::endl;
      } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "     eval = ((locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << "));"
              << std::endl;
          stream_program_src << "     index_calc = eval - (locIndex[(k*" << dims << ")+" << d << "]);" << std::endl;
          stream_program_src << "     abs = fabs(index_calc);" << std::endl;
          stream_program_src << "     last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
          stream_program_src << "     localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
              << std::endl;
          stream_program_src << "     curSupport *= localSupport;" << std::endl << std::endl;
        }
      }

      stream_program_src << "     myResult += curSupport;" << std::endl;
      stream_program_src << "  }" << std::endl;
      stream_program_src << std::endl;
      stream_program_src << "  barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream_program_src << " }" << std::endl;
      stream_program_src << std::endl;
      stream_program_src << " // Iterate over all grid points (slow ones, without cache)" << std::endl;
      stream_program_src << " for(int m = start_grid + fastChunkSizeGrid; m < end_grid; m++)" << std::endl;
      stream_program_src << " {" << std::endl;
      stream_program_src << "   curSupport = ptrAlpha[m];" << std::endl << std::endl;

      if (dims > maxDimUnroll) {
        stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "   eval = ((ptrLevel[(m*" << dims << ")+ d]) * (data[d]));" << std::endl;
        stream_program_src << "   index_calc = eval - (ptrIndex[(m*" << dims << ")+ d]);" << std::endl;
        stream_program_src << "   abs = fabs(index_calc);" << std::endl;
        stream_program_src << "   last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
        stream_program_src << "   localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
            << std::endl;
        stream_program_src << "   curSupport *= localSupport;" << std::endl << std::endl;
        stream_program_src << "}" << std::endl;
      } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "   eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));"
              << std::endl;
          stream_program_src << "   index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
          stream_program_src << "   abs = fabs(index_calc);" << std::endl;
          stream_program_src << "   last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
          stream_program_src << "   localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
              << std::endl;
          stream_program_src << "   curSupport *= localSupport;" << std::endl << std::endl;
        }
      }

      stream_program_src << "   myResult += curSupport;" << std::endl;
      stream_program_src << " }" << std::endl;
    } else {
      stream_program_src << " // Iterate over all grid points (without cache)" << std::endl;
      stream_program_src << " for(int m = start_grid; m < end_grid; m++)" << std::endl;
      stream_program_src << " {" << std::endl;
      stream_program_src << "   curSupport = ptrAlpha[m];" << std::endl << std::endl;

      if (dims > maxDimUnroll) {
        stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "   eval = ((ptrLevel[(m*" << dims << ") + d]) * (data[d]));" << std::endl;
        stream_program_src << "   index_calc = eval - (ptrIndex[(m*" << dims << ")+ d]);" << std::endl;
        stream_program_src << "   abs = fabs(index_calc);" << std::endl;
        stream_program_src << "   last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
        stream_program_src << "   localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
            << std::endl;
        stream_program_src << "   curSupport *= localSupport;" << std::endl << std::endl;
        //stream_program_src << "   printf(\"l: %lf, i: %lf, sur: %lf\\n\", ptrLevel[m*" << dims <<" + d], ptrIndex[m*" << dims <<" + d], ptrAlpha[m]);" << std::endl;

        stream_program_src << "}" << std::endl;
      } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "   eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));"
              << std::endl;
          stream_program_src << "   index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
          stream_program_src << "   abs = fabs(index_calc);" << std::endl;
          stream_program_src << "   last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
          stream_program_src << "   localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
              << std::endl;
          stream_program_src << "   curSupport *= localSupport;" << std::endl << std::endl;
        }
      }

      stream_program_src << "   myResult += curSupport;" << std::endl;
      stream_program_src << " }" << std::endl;

    }
    stream_program_src << std::endl;
    stream_program_src << " ptrResult[globalIdx] = myResult;" << std::endl;
    stream_program_src << "}" << std::endl;

    //std::cout << stream_program_src.str();

    //update file with kernel (for debugging)
    std::ofstream multFile;
    multFile.open("multKernel_tmp.cl");
    multFile << stream_program_src.str();
    multFile.close();

    return stream_program_src.str();
  }

  std::string generateSourceMultTrans(size_t dims) {

    size_t localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    bool useLocalMemory = this->parameters.getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY");
    uint64_t maxDimUnroll = this->parameters.getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL");

    std::stringstream stream_program_src;

    if (getType<real_type>::asString() == "double") {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    stream_program_src << "__kernel" << std::endl;

    stream_program_src << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    stream_program_src << "void multTransOCL(__global const " << getType<real_type>::asString() << "* ptrLevel,"
        << std::endl;
    stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrIndex," << std::endl;
    stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrData," << std::endl;
    stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrSource," << std::endl;
    stream_program_src << "           __global       " << getType<real_type>::asString() << "* ptrResult," << std::endl;
    stream_program_src << "           uint sourceSize," << std::endl;
    stream_program_src << "           uint start_data," << std::endl;
    stream_program_src << "           uint end_data)" << std::endl;
    stream_program_src << "{" << std::endl;
    stream_program_src << " int globalIdx = get_global_id(0);" << std::endl;
    stream_program_src << " int localIdx = get_local_id(0);" << std::endl;
    stream_program_src << std::endl;
    stream_program_src << " " << getType<real_type>::asString()
        << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
    stream_program_src << " " << getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl
        << std::endl;

    if (useLocalMemory) {
      stream_program_src << " __local " << getType<real_type>::asString() << " locData[" << dims * localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << " __local " << getType<real_type>::asString() << " locSource[" << localWorkgroupSize << "];"
          << std::endl << std::endl;
    }

    if (dims > maxDimUnroll) {
      stream_program_src << " " << getType<real_type>::asString() << " level[" << dims << "];" << std::endl;
      stream_program_src << " " << getType<real_type>::asString() << " index[" << dims << "];" << std::endl;
      stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
      stream_program_src << "  level[d] = ptrLevel[(globalIdx*" << dims << ") + d];" << std::endl;
      stream_program_src << "  index[d] = ptrIndex[(globalIdx*" << dims << ") + d];" << std::endl;
      stream_program_src << "}" << std::endl;
    } else {
      for (size_t d = 0; d < dims; d++) {
        stream_program_src << " " << getType<real_type>::asString() << " level_" << d << " = ptrLevel[(globalIdx*"
            << dims << ")+" << d << "];" << std::endl;
        stream_program_src << " " << getType<real_type>::asString() << " index_" << d << " = ptrIndex[(globalIdx*"
            << dims << ")+" << d << "];" << std::endl;
      }
    }

    stream_program_src << std::endl;
    stream_program_src << " // Iterate over all grid points" << std::endl;

    if (useLocalMemory) {
      stream_program_src << " for(int i = start_data; i < end_data; i+=" << localWorkgroupSize << ")" << std::endl;
      stream_program_src << " {" << std::endl;

      if (dims > maxDimUnroll) {
        stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "   locData[(d * " << localWorkgroupSize
            << ")+(localIdx)] = ptrData[(d * sourceSize) + (localIdx + i)];" << std::endl;
        stream_program_src << "}" << std::endl;
      } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "   locData[(" << d << "*" << localWorkgroupSize << ")+(localIdx)] = ptrData[(" << d
              << "*sourceSize)+(localIdx+i)];" << std::endl;
        }
      }

      stream_program_src << "   locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
      stream_program_src << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl << std::endl;
      stream_program_src << "   for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
      stream_program_src << "   {" << std::endl;

      stream_program_src << "     curSupport = locSource[k];" << std::endl << std::endl;
    } else {
      stream_program_src << "   for(int k = start_data; k < end_data; k++)" << std::endl;
      stream_program_src << "   {" << std::endl;
      stream_program_src << "     curSupport = ptrSource[k];" << std::endl << std::endl;
    }

    if (dims > maxDimUnroll) {
      stream_program_src << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;

      if (useLocalMemory) {
        stream_program_src << "     eval = (level[d] * (locData[(d * " << localWorkgroupSize << ")+k]));" << std::endl;
      } else {
        stream_program_src << "     eval = (level[d] * (ptrData[(d * sourceSize) + k]));" << std::endl;
      }

      stream_program_src << "     index_calc = eval - index[d];" << std::endl;
      stream_program_src << "     abs = fabs(index_calc);" << std::endl;
      stream_program_src << "     last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
      stream_program_src << "     localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
          << std::endl;
      stream_program_src << "     curSupport *= localSupport;" << std::endl;
      stream_program_src << "}" << std::endl;
    } else {
      for (size_t d = 0; d < dims; d++) {
        if (useLocalMemory) {
          stream_program_src << "     eval = ((level_" << d << ") * (locData[(" << d << "*" << localWorkgroupSize
              << ")+k]));" << std::endl;
        } else {
          stream_program_src << "     eval = ((level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]));"
              << std::endl;
        }

        stream_program_src << "     index_calc = eval - (index_" << d << ");" << std::endl;
        stream_program_src << "     abs = fabs(index_calc);" << std::endl;
        stream_program_src << "     last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
        stream_program_src << "     localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");"
            << std::endl;
        stream_program_src << "     curSupport *= localSupport;" << std::endl;
      }
    }

    stream_program_src << std::endl << "     myResult += curSupport;" << std::endl;
    stream_program_src << "   }" << std::endl << std::endl;
    if (useLocalMemory) {
      stream_program_src << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream_program_src << " }" << std::endl;
    }
    stream_program_src << " ptrResult[globalIdx] = myResult;" << std::endl;
    stream_program_src << "}" << std::endl;

//update file with kernel (for debugging)
    std::ofstream multTransFile;
    multTransFile.open("multTransKernel_tmp.cl");
    multTransFile << stream_program_src.str();
    multTransFile.close();

    return stream_program_src.str();
  }

}
;

}
}

