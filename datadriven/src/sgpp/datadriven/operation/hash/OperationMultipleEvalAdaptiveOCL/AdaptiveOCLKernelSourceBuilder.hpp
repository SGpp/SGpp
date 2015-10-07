/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>

//#include "AdaptiveOCLParameters.hpp"

namespace SGPP {
namespace datadriven {

namespace AdaptiveOCL {
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
class AdaptiveOCLKernelSourceBuilder {
private:
  std::shared_ptr<base::OCLConfigurationParameters> parameters;
  real_type blubb;
public:
  AdaptiveOCLKernelSourceBuilder(std::shared_ptr<base::OCLConfigurationParameters> parameters) :
      parameters(parameters) {

  }

  std::string generateSourceMult(size_t dims) {

    if (parameters->getAsBoolean("REUSE_SOURCE")) {
      std::stringstream stream_program_src;
      std::ifstream file;
      file.open("multKernel_tmp.cl");
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

    size_t localWorkgroupSize = parameters->getAsUnsigned("LOCAL_SIZE");
//    bool useLocalMemory = this->parameters->getAsBoolean("KERNEL_USE_LOCAL_MEMORY");
    uint64_t maxDimUnroll = this->parameters->getAsUnsigned("KERNEL_MAX_DIM_UNROLL");

    std::stringstream stream_program_src;

    std::string strType = AdaptiveOCL::getType<real_type>::asString();

    if (AdaptiveOCL::getType<real_type>::asString() == "double") {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    stream_program_src << "__kernel" << std::endl;
    stream_program_src << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    stream_program_src << "void multOCL(__global const " << strType << "* ptrStreamingArray, " << std::endl;
    stream_program_src << "           __global const " << strType << "* ptrSubspaceArray," << std::endl;
    stream_program_src << "           __global const uint* ptrMetaInfo," << std::endl;
    stream_program_src << "           __global const " << strType << "* ptrData," << std::endl;
    stream_program_src << "           __global const " << strType << "* ptrAlpha," << std::endl;
    stream_program_src << "           __global       " << strType << "* ptrResult," << std::endl;
    stream_program_src << "           uint resultSize," << std::endl;
    stream_program_src << "           uint numSubspaces)" << std::endl;
    stream_program_src << "{" << std::endl;
    stream_program_src << "  int globalIdx = get_global_id(0);" << std::endl;
    stream_program_src << "  int localIdx = get_local_id(0);" << std::endl;

    stream_program_src << std::endl;

    stream_program_src << "  " << strType << " index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
    stream_program_src << "  " << strType << " myResult = 0.0;" << std::endl << std::endl;
    stream_program_src << "  // Create registers for the data" << std::endl;

    if (dims > maxDimUnroll) {
      stream_program_src << "  " << strType << " data[" << dims << "];" << std::endl;
      stream_program_src << "  for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
      stream_program_src << "    data[d] = ptrData[globalIdx+(resultSize * d)];" << std::endl;
      stream_program_src << "  }" << std::endl;
    } else {
      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "  " << strType << " data_" << d
            << " = ptrData[globalIdx+(resultSize*" << d << ")];" << std::endl;
      }
    }

    stream_program_src << std::endl;

    stream_program_src << "  // Iterate over all subspaces (without cache)" << std::endl;
    stream_program_src << "  uint metaSize = "<< dims <<" + 3;" << std::endl;
    stream_program_src << "  for (uint l = 0; l < numSubspaces; l++)" << std::endl;
    stream_program_src << "  {" << std::endl;
    //create level variable
    if (dims > maxDimUnroll) {
        stream_program_src << "    " << strType << " level[" << dims << "];" << std::endl;
        stream_program_src << "    " << strType << " eval[" << dims << "];" << std::endl;
        stream_program_src << "    for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "      level[d] = (float)ptrMetaInfo[l*metaSize + d + 3];" << std::endl;
        stream_program_src << "      eval[d] = (pow(2,level[d]) * (data[d]));" << std::endl;
        stream_program_src << "    }" << std::endl;
     } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "    " << strType << " level_" << d
              << " = ptrMetaInfo[l*metaSize + " << d << " + 3];" << std::endl;
          stream_program_src << "    " << strType << " eval_" << d
                            << " = (pow(2.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ",level_" << d << ") * (data_" << d << "));" << std::endl;
     }
    }

    stream_program_src << "    uint type = ptrMetaInfo[l*metaSize];" << std::endl;
    stream_program_src << "    uint startIndex = ptrMetaInfo[l*metaSize + 1];" << std::endl;
    stream_program_src << "    uint numElements = ptrMetaInfo[l*metaSize + 2];" << std::endl;
    //STREAMING BRANCH
    stream_program_src << "    if ( type > 0 ) //stream" << std::endl;
    stream_program_src << "    {" << std::endl;
    stream_program_src << "        int streamElementSize = " << dims << " + 1;" << std::endl;
    stream_program_src << "        for ( int i = 0; i < numElements; i++)" << std::endl;
    stream_program_src << "        {" << std::endl;
    stream_program_src << "            int gridIndex = (int)ptrStreamingArray[streamElementSize*(startIndex + i) + " << dims <<"];" << std::endl;
    stream_program_src << "            curSupport = ptrAlpha[gridIndex];" << std::endl;
    stream_program_src << "            for (size_t d = 0; d < " << dims <<"; d++)" << std::endl;
    stream_program_src << "            {" << std::endl;
    stream_program_src << "                index_calc = eval[d] - ptrStreamingArray[streamElementSize*(startIndex + i) + d];" << std::endl;
    stream_program_src << "                abs = fabs(index_calc);" << std::endl;
    stream_program_src << "                last = 1.0" << AdaptiveOCL::getType<real_type>::constSuffix() << " - abs;" << std::endl;
    stream_program_src << "                localSupport = fmax(last, 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ");" << std::endl;
    stream_program_src << "                curSupport *= localSupport;" << std::endl;
    stream_program_src << "            }" << std::endl;
    stream_program_src << "            myResult += curSupport;" << std::endl;
    stream_program_src << "        }" << std::endl;
    stream_program_src << "    }" << std::endl;
    //SUBSPACE BRANCH
    stream_program_src << "    else //subspace" << std::endl;
    stream_program_src << "    {" << std::endl;
    stream_program_src << "        int subElementSize = " << dims << " + 1;" << std::endl;
    //RELEVANT INDEX
    stream_program_src << "        //calc relevant index" << std::endl;
    stream_program_src << "        int floored = 0;" << std::endl;
    stream_program_src << "        " << strType << " index[" << dims <<"];" << std::endl;
    stream_program_src << "        for (size_t d = 0; d < " << dims <<"; d++)" << std::endl;
    stream_program_src << "        {" << std::endl;
    stream_program_src << "            floored = (int)(data[d]*pow(2.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ",level[d]));" << std::endl;
    stream_program_src << "            index[d] = (" << strType << ")(floored | 0x1);" << std::endl;
    stream_program_src << "        }" << std::endl;
    //LINEAR INDEX
    stream_program_src << "        //calc linear index of gridpoint (relative to subspace)" << std::endl;
    stream_program_src << "        " << strType << " result = 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ";" << std::endl;
    stream_program_src << "        " << strType << " index_half = 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ";" << std::endl;
    stream_program_src << "        " << strType << " level_calc = 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ";" << std::endl;
    stream_program_src << "        for (size_t d = 0; d < "<< dims-1<< "; d++) {" << std::endl;
    stream_program_src << "            index_half = (int)index[d] >> 1;" << std::endl;
    stream_program_src << "            level_calc = pow(2.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ", level[d+1]-1);" << std::endl;
    stream_program_src << "            result += index_half;" << std::endl;
    stream_program_src << "            result *= level_calc;" << std::endl;
    stream_program_src << "        }" << std::endl;
    stream_program_src << "        result += (int)index["<< dims - 1 <<"] >> 1;" << std::endl;
    stream_program_src << "        size_t linIndex = (size_t)result;" << std::endl << std::endl;
    //GRID INDEX AND EVAL
    stream_program_src << "        " << strType << " gridIndex = ptrSubspaceArray[subElementSize*(startIndex + linIndex) + "<< dims <<"];" << std::endl;
    stream_program_src << "        if ( !isnan(gridIndex) )" << std::endl;
    stream_program_src << "        {" << std::endl;
    stream_program_src << "            curSupport = ptrAlpha[(int)gridIndex];" << std::endl << std::endl;

    stream_program_src << "            for (size_t d = 0; d < " << dims << "; d++)" << std::endl;
    stream_program_src << "            {" << std::endl;
    stream_program_src << "                index_calc = eval[d] - index[d];" << std::endl;
    stream_program_src << "                abs = fabs(index_calc);" << std::endl;
    stream_program_src << "                last = 1.0" << AdaptiveOCL::getType<real_type>::constSuffix() << " - abs;" << std::endl;
    stream_program_src << "                localSupport = fmax(last, 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ");" << std::endl;
    stream_program_src << "                curSupport *= localSupport;" << std::endl;
    stream_program_src << "            }" << std::endl;
    stream_program_src << "            myResult += curSupport;" << std::endl;
    stream_program_src << "        }" << std::endl;
    stream_program_src << "    }" << std::endl;

    stream_program_src << std::endl;
    stream_program_src << "  ptrResult[globalIdx] = myResult;" << std::endl;
    stream_program_src << "  }" << std::endl;
    stream_program_src << "}" << std::endl;

    //update file with kernel (for debugging)
    std::ofstream multFile;
    multFile.open("multKernel_tmp.cl");
    multFile << stream_program_src.str();
    multFile.close();

    return stream_program_src.str();
  }

  std::string generateSourceMultTrans(size_t dims) {

    if (parameters->getAsBoolean("REUSE_SOURCE")) {
      std::stringstream stream_program_src;
      std::ifstream file;
      file.open("multTransKernel_tmp.cl");
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

    size_t localWorkgroupSize = parameters->getAsUnsigned("LOCAL_SIZE");
//    bool useLocalMemory = this->parameters->getAsBoolean("KERNEL_USE_LOCAL_MEMORY");
    uint64_t maxDimUnroll = this->parameters->getAsUnsigned("KERNEL_MAX_DIM_UNROLL");

    std::stringstream stream_program_src;
    std::string strPrecisionType = AdaptiveOCL::getType<real_type>::asString();

    if (strPrecisionType == "double") {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl;
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable" << std::endl << std::endl;

      //atomic add for doubles
      stream_program_src << "void atomic_add_global(volatile __global double *source, const double operand) {" << std::endl;
      stream_program_src << "  union {" << std::endl;
      stream_program_src << "    long intVal;" << std::endl;
      stream_program_src << "    double floatVal;" << std::endl;
      stream_program_src << "  } newVal;" << std::endl;
      stream_program_src << "  union {" << std::endl;
      stream_program_src << "    long intVal;" << std::endl;
      stream_program_src << "    double floatVal;" << std::endl;
      stream_program_src << "  } prevVal;" << std::endl;
      stream_program_src << "  do {" << std::endl;
      stream_program_src << "    prevVal.floatVal = *source;" << std::endl;
      stream_program_src << "    newVal.floatVal = prevVal.floatVal + operand;" << std::endl;
      stream_program_src << "  } while (atom_cmpxchg((volatile __global long *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);" << std::endl;
      stream_program_src << "}" << std::endl << std::endl;
    }
    else
    {
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl;
      stream_program_src << "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable" << std::endl << std::endl;

      //atomic add for float
      stream_program_src << "void atomic_add_global(volatile __global float *source, const float operand) {" << std::endl;
      stream_program_src << "  union {" << std::endl;
      stream_program_src << "    uint intVal;" << std::endl;
      stream_program_src << "    float floatVal;" << std::endl;
      stream_program_src << "  } newVal;" << std::endl;
      stream_program_src << "  union {" << std::endl;
      stream_program_src << "    uint intVal;" << std::endl;
      stream_program_src << "    float floatVal;" << std::endl;
      stream_program_src << "  } prevVal;" << std::endl;
      stream_program_src << "  do {" << std::endl;
      stream_program_src << "    prevVal.floatVal = *source;" << std::endl;
      stream_program_src << "    newVal.floatVal = prevVal.floatVal + operand;" << std::endl;
      stream_program_src << "  } while (atomic_cmpxchg((volatile __global uint *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);" << std::endl;
      stream_program_src << "}" << std::endl << std::endl;
    }

    stream_program_src << "__kernel" << std::endl;

    stream_program_src << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    stream_program_src << "void multTransOCL(__global const " << strPrecisionType << "* ptrStreamingArray," << std::endl;
    stream_program_src << "           __global const " << strPrecisionType << "* ptrSubspaceArray," << std::endl;
    stream_program_src << "           __global const uint* ptrMetaInfo," << std::endl;
    stream_program_src << "           __global const " << strPrecisionType << "* ptrData," << std::endl;
    stream_program_src << "           __global const " << strPrecisionType << "* ptrSource," << std::endl;
    stream_program_src << "           __global       " << strPrecisionType << "* ptrResult," << std::endl;
    stream_program_src << "           uint sourceSize," << std::endl;
    stream_program_src << "           uint numSubspaces)" << std::endl;
    stream_program_src << "{" << std::endl;
    stream_program_src << "    int globalIdx = get_global_id(0);" << std::endl;
    stream_program_src << "    int localIdx = get_local_id(0);" << std::endl;
    stream_program_src << std::endl;
    stream_program_src << "    " << strPrecisionType << " index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;


    if (dims > maxDimUnroll) {
      stream_program_src << "    " << strPrecisionType << " data[" << dims << "];" << std::endl;
      stream_program_src << "    for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
      stream_program_src << "        data[d] = ptrData[globalIdx+(sourceSize * d)];" << std::endl;
      stream_program_src << "    }" << std::endl;
    } else {
      for (size_t d = 0; d < dims; d++) {
        stream_program_src << "    " << strPrecisionType << " data_" << d << " = ptrData[globalIdx+(sourceSize *"<< d << ")];" << std::endl;
      }
    }

    stream_program_src << std::endl;
    stream_program_src << "    // Iterate over all subspaces" << std::endl;
    stream_program_src << "    uint metaSize = "<< dims <<" + 3;" << std::endl;
    stream_program_src << "    for (uint l = 0; l < numSubspaces; l++)" << std::endl;
    stream_program_src << "    {" << std::endl;
    //create level variable
    if (dims > maxDimUnroll) {
        stream_program_src << "        " << strPrecisionType << " level[" << dims << "];" << std::endl;
        stream_program_src << "        " << strPrecisionType << " eval[" << dims << "];" << std::endl;
        stream_program_src << "        for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        stream_program_src << "            level[d] = (float)ptrMetaInfo[l*metaSize + d + 3];" << std::endl;
        stream_program_src << "            eval[d] = (pow(2,level[d]) * (data[d]));" << std::endl;
        stream_program_src << "        }" << std::endl;
     } else {
        for (size_t d = 0; d < dims; d++) {
          stream_program_src << "        " << strPrecisionType << " level_" << d
              << " = ptrMetaInfo[l*metaSize + " << d << " + 3];" << std::endl;
          stream_program_src << "    " << strPrecisionType << " eval_" << d
                            << " = (pow(2,level_" << d << ") * (data_" << d << "));" << std::endl;
     }
    }


    stream_program_src << "        uint type = ptrMetaInfo[l*metaSize];" << std::endl;
    stream_program_src << "        uint startIndex = ptrMetaInfo[l*metaSize + 1];" << std::endl;
    stream_program_src << "        uint numElements = ptrMetaInfo[l*metaSize + 2];" << std::endl;
    //STREAMING BRANCH
    stream_program_src << "        if ( type > 0 ) //stream" << std::endl;
    stream_program_src << "        {" << std::endl;
    stream_program_src << "            int streamElementSize = " << dims << " + 1;" << std::endl;
    stream_program_src << "            for ( int i = 0; i < numElements; i++)" << std::endl;
    stream_program_src << "            {" << std::endl;

    stream_program_src << "                curSupport = ptrSource[globalIdx];" << std::endl;
    stream_program_src << "                for (size_t d = 0; d < " << dims <<"; d++)" << std::endl;
    stream_program_src << "                {" << std::endl;
    stream_program_src << "                    index_calc = eval[d] - ptrStreamingArray[streamElementSize*(startIndex + i) + d];" << std::endl;
    stream_program_src << "                    abs = fabs(index_calc);" << std::endl;
    stream_program_src << "                    last = 1.0" << AdaptiveOCL::getType<real_type>::constSuffix() << " - abs;" << std::endl;
    stream_program_src << "                    localSupport = fmax(last, 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ");" << std::endl;
    stream_program_src << "                    curSupport *= localSupport;" << std::endl;
    stream_program_src << "                }" << std::endl;
    stream_program_src << "                int gridIndex = (int)ptrStreamingArray[streamElementSize*(startIndex + i) + " << dims <<"];" << std::endl;
    stream_program_src << "                atomic_add_global(&ptrResult[gridIndex], curSupport);" << std::endl;
    stream_program_src << "            }" << std::endl;
    stream_program_src << "        }" << std::endl;
    //SUBSPACE BRANCH
    stream_program_src << "        else //subspace" << std::endl;
    stream_program_src << "        {" << std::endl;
    stream_program_src << "            int subElementSize = " << dims << " + 1;" << std::endl;
    //RELEVANT INDEX
    stream_program_src << "            //calc relevant index" << std::endl;
    stream_program_src << "            int floored = 0;" << std::endl;
    stream_program_src << "            " << strPrecisionType << " index[" << dims <<"];" << std::endl;
    stream_program_src << "            for (size_t d = 0; d < " << dims <<"; d++)" << std::endl;
    stream_program_src << "            {" << std::endl;
    stream_program_src << "                floored = (int)(data[d]*pow(2.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ",level[d]));" << std::endl;
    stream_program_src << "                index[d] = (" << strPrecisionType << ")(floored | 0x1);" << std::endl;
    stream_program_src << "            }" << std::endl;
    //LINEAR INDEX
    stream_program_src << "            //calc linear index of gridpoint (relative to subspace)" << std::endl;
    stream_program_src << "            " << strPrecisionType << " result = 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ";" << std::endl;
    stream_program_src << "            " << strPrecisionType << " index_half = 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ";" << std::endl;
    stream_program_src << "            " << strPrecisionType << " level_calc = 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ";" << std::endl;
    stream_program_src << "            for (size_t d = 0; d < "<< dims-1<< "; d++) {" << std::endl;
    stream_program_src << "                index_half = (int)index[d] >> 1;" << std::endl;
    stream_program_src << "                level_calc = pow(2.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ", level[d+1]-1);" << std::endl;
    stream_program_src << "                result += index_half;" << std::endl;
    stream_program_src << "                result *= level_calc;" << std::endl;
    stream_program_src << "            }" << std::endl;
    stream_program_src << "            result += (int)index["<< dims - 1 <<"] >> 1;" << std::endl;
    stream_program_src << "            size_t linIndex = (size_t)result;" << std::endl << std::endl;
    //GRID INDEX AND EVAL
    stream_program_src << "            "<< strPrecisionType <<" gridIndex = ptrSubspaceArray[subElementSize*(startIndex + linIndex) + " << dims << "];" << std::endl;
    stream_program_src << "            if ( !isnan(gridIndex) )" << std::endl;
    stream_program_src << "            {" << std::endl;
    stream_program_src << "                curSupport = ptrSource[globalIdx];" << std::endl << std::endl;
    stream_program_src << "                for (size_t d = 0; d < " << dims << "; d++)" << std::endl;
    stream_program_src << "                {" << std::endl;
    stream_program_src << "                    index_calc = eval[d] - index[d];" << std::endl;
    stream_program_src << "                    abs = fabs(index_calc);" << std::endl;
    stream_program_src << "                    last = 1.0" << AdaptiveOCL::getType<real_type>::constSuffix() << " - abs;" << std::endl;
    stream_program_src << "                    localSupport = fmax(last, 0.0" << AdaptiveOCL::getType<real_type>::constSuffix() << ");" << std::endl;
    stream_program_src << "                    curSupport *= localSupport;" << std::endl;
    stream_program_src << "                }" << std::endl;
    stream_program_src << "                atomic_add_global(&ptrResult[(int)gridIndex], curSupport);" << std::endl;
    stream_program_src << "            }" << std::endl;
    stream_program_src << "        }" << std::endl;
    stream_program_src << "    }" << std::endl;
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

