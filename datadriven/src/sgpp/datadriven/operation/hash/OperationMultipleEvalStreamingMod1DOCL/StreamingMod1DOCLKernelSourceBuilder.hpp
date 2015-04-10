/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>

//#include "StreamingOCLParameters.hpp"

namespace SGPP {
namespace datadriven {

namespace StreamingMod1DOCL {
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
}

template<typename real_type>
class StreamingMod1DOCLKernelSourceBuilder {
private:
  base::ConfigurationParameters parameters;
  real_type blubb;
public:
  StreamingMod1DOCLKernelSourceBuilder(base::ConfigurationParameters parameters) :
      parameters(parameters) {

  }

  std::string generateSourceMult(size_t dims) {

    if (parameters.getAsBoolean("REUSE_SOURCE")) {
      std::stringstream stream_program_src;
      std::ifstream file;
      file.open("multKernelMod1D_tmp.cl");
      if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {
          stream_program_src << line << std::endl;
        }
        file.close();
      } else {
        throw new base::operation_exception("OCL error: file to reuse not found\n");
      }
      return stream_program_src.str();
    }

    size_t localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    bool useLocalMemory = this->parameters.getAsBoolean("KERNEL_USE_LOCAL_MEMORY");

    std::stringstream stream_program_src;

    if (StreamingMod1DOCL::getType<real_type>::asString() == "double") {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    stream_program_src << "__kernel" << std::endl;
    stream_program_src << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    stream_program_src << "void multOCL(__global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrLevel,"
        << std::endl;
    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrIndex," << std::endl;
//    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrMask," << std::endl; // not needed for this kernel, but there for uniformity
//    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrOffset," << std::endl; // not needed for this kernel, but there for uniformity
    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrData," << std::endl;
    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrAlpha," << std::endl;
    stream_program_src << "           __global       " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrResult," << std::endl;
    stream_program_src << "           uint resultSize," << std::endl;
    stream_program_src << "           uint start_grid," << std::endl;
    stream_program_src << "           uint end_grid) " << std::endl;
    stream_program_src << "{" << std::endl;
    stream_program_src << "   int globalIdx = get_global_id(0);" << std::endl;
    stream_program_src << "   int localIdx = get_local_id(0);" << std::endl;
    stream_program_src << std::endl;
    if (useLocalMemory) {
      stream_program_src << "   __local " << StreamingMod1DOCL::getType<real_type>::asString() << " locLevel[" << dims * localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << "   __local " << StreamingMod1DOCL::getType<real_type>::asString() << " locIndex[" << dims * localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << "   __local " << StreamingMod1DOCL::getType<real_type>::asString() << " locAlpha[" << localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << std::endl;
    }
    stream_program_src << "   " << StreamingMod1DOCL::getType<real_type>::asString()
        << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
    stream_program_src << "   " << StreamingMod1DOCL::getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl
        << std::endl;
    stream_program_src << "   // Create registers for the data" << std::endl;

    for (size_t d = 0; d < dims; d++) {
      stream_program_src << " " << StreamingMod1DOCL::getType<real_type>::asString() << " data_" << d
          << " = ptrData[globalIdx+(resultSize*" << d << ")];" << std::endl;
    }

    stream_program_src << std::endl;
    if (useLocalMemory) {
      stream_program_src << "   // Iterate over all grid points (fast ones, with cache)" << std::endl;
      stream_program_src << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
      stream_program_src << " uint fastChunkSizeGrid = (chunkSizeGrid / " << localWorkgroupSize << ") * "
          << localWorkgroupSize << ";" << std::endl;
      stream_program_src << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+=" << localWorkgroupSize
          << ")" << std::endl;
      stream_program_src << "   {" << std::endl;

      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "     locLevel[(localIdx*" << dims << ")+" << d << "] = ptrLevel[((j+localIdx)*" << dims
            << ")+" << d << "];" << std::endl;
        stream_program_src << "     locIndex[(localIdx*" << dims << ")+" << d << "] = ptrIndex[((j+localIdx)*" << dims
            << ")+" << d << "];" << std::endl;
      }

      stream_program_src << "       locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
      stream_program_src << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream_program_src << std::endl;
      stream_program_src << "       for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
      stream_program_src << "       {" << std::endl;
      stream_program_src << "           curSupport = locAlpha[k];" << std::endl << std::endl;

      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "         if ((locLevel[(k*" << dims << ")+" << d << "]) == 2.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")" << std::endl;
        stream_program_src << "         {" << std::endl;
        stream_program_src << "             curSupport *= 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ";" << std::endl;
        stream_program_src << "         }" << std::endl;
        stream_program_src << "         else if ((locIndex[(k*" << dims << ")+" << d << "]) == 1.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")" << std::endl;
        stream_program_src << "         {" << std::endl;
        stream_program_src << "             curSupport *= max(2.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - ( (locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << ") ), 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") ;" << std::endl;
        stream_program_src << "         }" << std::endl;
        stream_program_src << "         else if ((locIndex[(k*" << dims << ")+" << d << "]) == ((locLevel[(k*" << dims
            << ")+" << d << "]) - 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") )" << std::endl;
        stream_program_src << "         {" << std::endl;
        stream_program_src << "             curSupport *= max(( (locLevel[(k*" << dims << ")+" << d << "]) * (data_"
            << d << ") ) - (locIndex[(k*" << dims << ")+" << d << "]) + 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << ", 0.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
        stream_program_src << "         }" << std::endl;
        stream_program_src << "         else " << std::endl;
        stream_program_src << "         {" << std::endl;
        stream_program_src << "             curSupport *= max(1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - fabs( ( (locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (locIndex[(k*" << dims
            << ")+" << d << "]) ), 0.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
        stream_program_src << "         }" << std::endl;
      }

      stream_program_src << "           myResult += curSupport;" << std::endl;
      stream_program_src << "       }" << std::endl;
      stream_program_src << std::endl;
      stream_program_src << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream_program_src << "   }" << std::endl;
      stream_program_src << std::endl;
      stream_program_src << "   // Iterate over all grid points (slow ones, without cache)" << std::endl;
      stream_program_src << " for(int m = start_grid + fastChunkSizeGrid; m < end_grid; m++)" << std::endl;
      stream_program_src << "   {" << std::endl;
      stream_program_src << "       curSupport = ptrAlpha[m];" << std::endl << std::endl;

      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "     if ((ptrLevel[(m*" << dims << ")+" << d << "]) == 2.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")" << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ";" << std::endl;
        stream_program_src << "     }" << std::endl;
        stream_program_src << "     else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == 1.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")" << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= max(2.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ), 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") ;" << std::endl;
        stream_program_src << "     }" << std::endl;
        stream_program_src << "     else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == ((ptrLevel[(m*" << dims
            << ")+" << d << "]) - 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") )" << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= max(( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d
            << ") ) - (ptrIndex[(m*" << dims << ")+" << d << "]) + 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ", 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
        stream_program_src << "     }" << std::endl;
        stream_program_src << "     else " << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= max(1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - fabs( ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (ptrIndex[(m*" << dims
            << ")+" << d << "]) ), 0.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
        stream_program_src << "     }" << std::endl;
      }

      stream_program_src << "       myResult += curSupport;" << std::endl;
      stream_program_src << "   }" << std::endl;
    } else {
      stream_program_src << "   // Iterate over all grid points (slow ones, without cache)" << std::endl;
      stream_program_src << " for(int m = start_grid; m < end_grid; m++)" << std::endl;
      stream_program_src << "   {" << std::endl;
      stream_program_src << "       curSupport = ptrAlpha[m];" << std::endl << std::endl;

      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "     if ((ptrLevel[(m*" << dims << ")+" << d << "]) == 2.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")" << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ";" << std::endl;
        stream_program_src << "     }" << std::endl;
        stream_program_src << "     else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == 1.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")" << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= max(2.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ), 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") ;" << std::endl;
        stream_program_src << "     }" << std::endl;
        stream_program_src << "     else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == ((ptrLevel[(m*" << dims
            << ")+" << d << "]) - 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") )" << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= max(( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d
            << ") ) - (ptrIndex[(m*" << dims << ")+" << d << "]) + 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ", 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
        stream_program_src << "     }" << std::endl;
        stream_program_src << "     else " << std::endl;
        stream_program_src << "     {" << std::endl;
        stream_program_src << "         curSupport *= max(1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - fabs( ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (ptrIndex[(m*" << dims
            << ")+" << d << "]) ), 0.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
        stream_program_src << "     }" << std::endl;
      }

      stream_program_src << "       myResult += curSupport;" << std::endl;
      stream_program_src << "   }" << std::endl;
    }
    stream_program_src << std::endl;
    stream_program_src << "   ptrResult[globalIdx] = myResult;" << std::endl;
    stream_program_src << "}" << std::endl;

    //update file with kernel (for debugging)
    std::ofstream multFile;
    multFile.open("multKernelMod1D_tmp.cl");
    multFile << stream_program_src.str();
    multFile.close();

    return stream_program_src.str();
  }

  std::string generateSourceMultTrans(size_t dims) {

    if (parameters.getAsBoolean("REUSE_SOURCE")) {
      std::stringstream stream_program_src;
      std::ifstream file;
      file.open("multTransKernelMod1D_tmp.cl");
      if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {
          stream_program_src << line << std::endl;
        }
        file.close();
      } else {
        throw new base::operation_exception("OCL error: file to reuse not found\n");
      }
      return stream_program_src.str();
    }

    size_t localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    bool useLocalMemory = this->parameters.getAsBoolean("KERNEL_USE_LOCAL_MEMORY");

    std::stringstream stream_program_src;

    if (StreamingMod1DOCL::getType<real_type>::asString() == "double") {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    stream_program_src << "__kernel" << std::endl;
    stream_program_src << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    stream_program_src << "void multTransOCL(__global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrLevel,"
        << std::endl;
    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrIndex," << std::endl;
//    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrMask," << std::endl; // not needed for this kernel, but there for uniformity
//    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrOffset," << std::endl; // not needed for this kernel, but there for uniformity
    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrData," << std::endl;
    stream_program_src << "           __global const " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrSource," << std::endl;
    stream_program_src << "           __global       " << StreamingMod1DOCL::getType<real_type>::asString() << "* ptrResult," << std::endl;
    stream_program_src << "           uint sourceSize," << std::endl;
    stream_program_src << "           uint start_data," << std::endl;
    stream_program_src << "           uint end_data)" << std::endl;
    stream_program_src << "{" << std::endl;
    stream_program_src << "   int globalIdx = get_global_id(0);" << std::endl;
    stream_program_src << "   int localIdx = get_local_id(0);" << std::endl;
    stream_program_src << std::endl;
    stream_program_src << "   " << StreamingMod1DOCL::getType<real_type>::asString()
        << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
    stream_program_src << "   " << StreamingMod1DOCL::getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl
        << std::endl;
    if (useLocalMemory) {
      stream_program_src << "   __local " << StreamingMod1DOCL::getType<real_type>::asString() << " locData[" << dims * localWorkgroupSize
          << "];" << std::endl;
      stream_program_src << "   __local " << StreamingMod1DOCL::getType<real_type>::asString() << " locSource[" << localWorkgroupSize
          << "];" << std::endl << std::endl;
    }

    for (size_t d = 0; d < dims; d++) {
      stream_program_src << " " << StreamingMod1DOCL::getType<real_type>::asString() << " level_" << d << " = ptrLevel[(globalIdx*" << dims
          << ")+" << d << "];" << std::endl;
      stream_program_src << " " << StreamingMod1DOCL::getType<real_type>::asString() << " index_" << d << " = ptrIndex[(globalIdx*" << dims
          << ")+" << d << "];" << std::endl;
    }

    stream_program_src << std::endl;
    if (useLocalMemory) {
      stream_program_src << "   // Iterate over all grid points" << std::endl;
      stream_program_src << " for(int i = start_data; i < end_data; i+=" << localWorkgroupSize << ")" << std::endl;
      stream_program_src << "   {" << std::endl;

      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "     locData[(" << d << "*" << localWorkgroupSize << ")+(localIdx)] = ptrData[(" << d
            << "*sourceSize)+(localIdx+i)];" << std::endl;
      }

      stream_program_src << "       locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
      stream_program_src << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl << std::endl;
      stream_program_src << "       for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
      stream_program_src << "       {" << std::endl;
      stream_program_src << "           curSupport = locSource[k];" << std::endl << std::endl;
    } else {
      stream_program_src << "   for(int k = start_data; k < end_data; k++)" << std::endl;
      stream_program_src << "       {" << std::endl;
      stream_program_src << "           curSupport = ptrSource[k];" << std::endl << std::endl;
    }

    for (size_t d = 0; d < dims; d++) {
      stream_program_src << "         if ((level_" << d << ") == 2.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")"
          << std::endl;
      stream_program_src << "         {" << std::endl;
      stream_program_src << "             curSupport *= 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ";" << std::endl;
      stream_program_src << "         }" << std::endl;
      stream_program_src << "         else if ((index_" << d << ") == 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ")"
          << std::endl;
      stream_program_src << "         {" << std::endl;
      if (useLocalMemory) {
        stream_program_src << "             curSupport *= max(2.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - ( (level_" << d << ") * (locData[(" << d << "*" << localWorkgroupSize << ")+k]) ), 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") ;" << std::endl;
      } else {
        stream_program_src << "             curSupport *= max(2.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - ( (level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]) ), 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") ;" << std::endl;
      }
      stream_program_src << "         }" << std::endl;
      stream_program_src << "         else if ((index_" << d << ") == ((level_" << d << ") - 1.0"
          << StreamingMod1DOCL::getType<real_type>::constSuffix() << ") )" << std::endl;
      stream_program_src << "         {" << std::endl;
      if (useLocalMemory) {
        stream_program_src << "             curSupport *= max(( (level_" << d << ") * (locData[(" << d << "*"
            << localWorkgroupSize << ")+k]) ) - (index_" << d << ") + 1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << ", 0.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
      } else {
        stream_program_src << "             curSupport *= max(( (level_" << d << ") * (ptrData[(" << d
            << "*sourceSize)+k]) ) - (index_" << d << ") + 1.0, 0.0);" << std::endl;
      }
      stream_program_src << "         }" << std::endl;
      stream_program_src << "         else " << std::endl;
      stream_program_src << "         {" << std::endl;
      if (useLocalMemory) {
        stream_program_src << "             curSupport *= max(1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - fabs( ( (level_" << d << ") * (locData[(" << d << "*" << localWorkgroupSize << ")+k]) ) - (index_"
            << d << ") ), 0.0" << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
      } else {
        stream_program_src << "             curSupport *= max(1.0" << StreamingMod1DOCL::getType<real_type>::constSuffix()
            << " - fabs( ( (level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]) ) - (index_" << d << ") ), 0.0"
            << StreamingMod1DOCL::getType<real_type>::constSuffix() << ");" << std::endl;
      }
      stream_program_src << "         }" << std::endl;
    }

    stream_program_src << std::endl << "      myResult += curSupport;" << std::endl;
    stream_program_src << "       }" << std::endl << std::endl;
    if (useLocalMemory) {
      stream_program_src << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream_program_src << "   }" << std::endl;
    }
    stream_program_src << "   ptrResult[globalIdx] = myResult;" << std::endl;
    stream_program_src << "}" << std::endl;

    //update file with kernel (for debugging)
    std::ofstream multTransFile;
    multTransFile.open("multTransKernelMod1D_tmp.cl");
    multTransFile << stream_program_src.str();
    multTransFile.close();

    return stream_program_src.str();
  }

};

}
}

