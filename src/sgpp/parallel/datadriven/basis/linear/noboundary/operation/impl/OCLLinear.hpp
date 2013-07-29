/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OCLLINEAR_HPP
#define OCLLINEAR_HPP

#include <sstream>
#include "parallel/datadriven/basis/common/ocl/OCLKernelBase.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace parallel {
    template<typename real_type>
    class OCLLinear : public OCLKernelBase {
      public:
        static const KernelType kernelType = Standard;
      private:
        virtual std::string generateSourceMult(size_t dims, size_t local_workgroup_size) {
          std::stringstream stream_program_src;

          if (getType<real_type>::asString() == "double") {
            stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
          }

          stream_program_src << "__kernel" << std::endl;
          stream_program_src << "__attribute__((reqd_work_group_size(" << local_workgroup_size << ", 1, 1)))" << std::endl;
          stream_program_src << "void multOCL(__global const " << getType<real_type>::asString() << "* ptrLevel," << std::endl;
          stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrIndex," << std::endl;
          stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrMask," << std::endl; // not needed for this kernel, but there for uniformity
          stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrOffset," << std::endl; // not needed for this kernel, but there for uniformity
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
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << " __local " << getType<real_type>::asString() << " locLevel[" << dims* local_workgroup_size << "];" << std::endl;
          stream_program_src << " __local " << getType<real_type>::asString() << " locIndex[" << dims* local_workgroup_size << "];" << std::endl;
          stream_program_src << " __local " << getType<real_type>::asString() << " locAlpha[" << local_workgroup_size << "];" << std::endl;
          stream_program_src << std::endl;
#endif
          stream_program_src << " " << getType<real_type>::asString() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
          stream_program_src << " " << getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl << std::endl;
          stream_program_src << " // Create registers for the data" << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << " " << getType<real_type>::asString() << " data_" << d << " = ptrData[globalIdx+(resultSize*" << d << ")];" << std::endl;
          }

          stream_program_src << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << " // Iterate over all grid points (fast ones, with cache)" << std::endl;
          stream_program_src << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
          stream_program_src << " uint fastChunkSizeGrid = (chunkSizeGrid / " << local_workgroup_size << ") * " << local_workgroup_size << ";" << std::endl;
          stream_program_src << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+=" << local_workgroup_size << ")" << std::endl;
          stream_program_src << " {" << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "   locLevel[(localIdx*" << dims << ")+" << d << "] = ptrLevel[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
            stream_program_src << "   locIndex[(localIdx*" << dims << ")+" << d << "] = ptrIndex[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
          }

          stream_program_src << "   locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
          stream_program_src << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          stream_program_src << std::endl;
          stream_program_src << "   for(int k = 0; k < " << local_workgroup_size << "; k++)" << std::endl;
          stream_program_src << "   {" << std::endl;
          stream_program_src << "     curSupport = locAlpha[k];" << std::endl << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "     eval = ((locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << "));" << std::endl;
            stream_program_src << "     index_calc = eval - (locIndex[(k*" << dims << ")+" << d << "]);" << std::endl;
            stream_program_src << "     abs = fabs(index_calc);" << std::endl;
            stream_program_src << "     last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
            stream_program_src << "     localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "     curSupport *= localSupport;" << std::endl << std::endl;
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

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "   eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));" << std::endl;
            stream_program_src << "   index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
            stream_program_src << "   abs = fabs(index_calc);" << std::endl;
            stream_program_src << "   last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
            stream_program_src << "   localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "   curSupport *= localSupport;" << std::endl << std::endl;
          }

          stream_program_src << "   myResult += curSupport;" << std::endl;
          stream_program_src << " }" << std::endl;
#else
          stream_program_src << " // Iterate over all grid points (without cache)" << std::endl;
          stream_program_src << " for(int m = start_grid; m < end_grid; m++)" << std::endl;
          stream_program_src << " {" << std::endl;
          stream_program_src << "   curSupport = ptrAlpha[m];" << std::endl << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "   eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));" << std::endl;
            stream_program_src << "   index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
            stream_program_src << "   abs = fabs(index_calc);" << std::endl;
            stream_program_src << "   last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
            stream_program_src << "   localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "   curSupport *= localSupport;" << std::endl << std::endl;
          }

          stream_program_src << "   myResult += curSupport;" << std::endl;
          stream_program_src << " }" << std::endl;

#endif
          stream_program_src << std::endl;
          stream_program_src << " ptrResult[globalIdx] = myResult;" << std::endl;
          stream_program_src << "}" << std::endl;

          return stream_program_src.str();
        }

        virtual std::string generateSourceMultTrans(size_t dims, size_t local_workgroup_size) {
          std::stringstream stream_program_src;

          if (getType<real_type>::asString() == "double") {
            stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
          }

          stream_program_src << "__kernel" << std::endl;;
          stream_program_src << "__attribute__((reqd_work_group_size(" << local_workgroup_size << ", 1, 1)))" << std::endl;
          stream_program_src << "void multTransOCL(__global const " << getType<real_type>::asString() << "* ptrLevel," << std::endl;
          stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrIndex," << std::endl;
          stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrMask," << std::endl; // not needed for this kernel, but there for uniformity
          stream_program_src << "           __global const " << getType<real_type>::asString() << "* ptrOffset," << std::endl; // not needed for this kernel, but there for uniformity
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
          stream_program_src << " " << getType<real_type>::asString() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
          stream_program_src << " " << getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << " __local " << getType<real_type>::asString() << " locData[" << dims* local_workgroup_size << "];" << std::endl;
          stream_program_src << " __local " << getType<real_type>::asString() << " locSource[" << local_workgroup_size << "];" << std::endl << std::endl;
#endif

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << " " << getType<real_type>::asString() << " level_" << d << " = ptrLevel[(globalIdx*" << dims << ")+" << d << "];" << std::endl;
            stream_program_src << " " << getType<real_type>::asString() << " index_" << d << " = ptrIndex[(globalIdx*" << dims << ")+" << d << "];" << std::endl;
          }

          stream_program_src << std::endl;
          stream_program_src << " // Iterate over all grid points" << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << " for(int i = start_data; i < end_data; i+=" << local_workgroup_size << ")" << std::endl;
          stream_program_src << " {" << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "   locData[(" << d << "*" << local_workgroup_size << ")+(localIdx)] = ptrData[(" << d << "*sourceSize)+(localIdx+i)];" << std::endl;
          }

          stream_program_src << "   locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
          stream_program_src << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl << std::endl;
          stream_program_src << "   for(int k = 0; k < " << local_workgroup_size << "; k++)" << std::endl;
          stream_program_src << "   {" << std::endl;

          stream_program_src << "     curSupport = locSource[k];" << std::endl << std::endl;
#else
          stream_program_src << "   for(int k = start_data; k < end_data; k++)" << std::endl;
          stream_program_src << "   {" << std::endl;
          stream_program_src << "     curSupport = ptrSource[k];" << std::endl << std::endl;

#endif

          for (size_t d = 0; d < dims; d++) {
#ifdef USEOCL_LOCAL_MEMORY
            stream_program_src << "     eval = ((level_" << d << ") * (locData[(" << d << "*" << local_workgroup_size << ")+k]));" << std::endl;
#else
            stream_program_src << "     eval = ((level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]));" << std::endl;
#endif
            stream_program_src << "     index_calc = eval - (index_" << d << ");" << std::endl;
            stream_program_src << "     abs = fabs(index_calc);" << std::endl;
            stream_program_src << "     last = 1.0" << getType<real_type>::constSuffix() << " - abs;" << std::endl;
            stream_program_src << "     localSupport = fmax(last, 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "     curSupport *= localSupport;" << std::endl;
          }

          stream_program_src << std::endl << "     myResult += curSupport;" << std::endl;
          stream_program_src << "   }" << std::endl << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          stream_program_src << " }" << std::endl;
#endif
          stream_program_src << " ptrResult[globalIdx] = myResult;" << std::endl;
          stream_program_src << "}" << std::endl;

          return stream_program_src.str();
        }
      public:
        static inline void multDefault(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* /*mask*/, //unused for this specialization
          sg::base::DataMatrix* /*offset*/, //unused for this specialization
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& alpha,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
        static inline void multDefault(
          sg::base::DataMatrixSP* level,
          sg::base::DataMatrixSP* index,
          sg::base::DataMatrixSP* /*mask*/, //unused for this specialization
          sg::base::DataMatrixSP* /*offset*/, //unused for this specialization
          sg::base::DataMatrixSP* dataset,
          sg::base::DataVectorSP& alpha,
          sg::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
        static inline void multTransposeDefault(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* /*mask*/, //unused for this specialization
          sg::base::DataMatrix* /*offset*/, //unused for this specialization
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& source,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
        static inline void multTransposeDefault(
          sg::base::DataMatrixSP* level,
          sg::base::DataMatrixSP* index,
          sg::base::DataMatrixSP* /*mask*/, //unused for this specialization
          sg::base::DataMatrixSP* /*offset*/, //unused for this specialization
          sg::base::DataMatrixSP* dataset,
          sg::base::DataVectorSP& source,
          sg::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
    };

    template<>
    inline void OCLLinear<double>::multDefault(
      sg::base::DataMatrix* level,
      sg::base::DataMatrix* index,
      sg::base::DataMatrix* /*mask*/, //unused for this specialization
      sg::base::DataMatrix* /*offset*/, //unused for this specialization
      sg::base::DataMatrix* dataset,
      sg::base::DataVector& alpha,
      sg::base::DataVector& result,
      const size_t start_index_grid,
      const size_t end_index_grid,
      const size_t start_index_data,
      const size_t end_index_data) {
      double* ptrLevel = level->getPointer();
      double* ptrIndex = index->getPointer();
      double* ptrAlpha = alpha.getPointer();
      double* ptrData = dataset->getPointer();
      double* ptrResult = result.getPointer();
      size_t result_size = result.getSize();
      size_t dims = dataset->getNrows();

      for (size_t i = start_index_data; i < end_index_data; i++) {
        for (size_t j = start_index_grid; j < end_index_grid; j++) {
          double curSupport = ptrAlpha[j];

          for (size_t d = 0; d < dims; d++) {
            double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
            double index_calc = eval - (ptrIndex[(j * dims) + d]);
            double abs = fabs(index_calc);
            double last = 1.0 - abs;
            double localSupport = std::max<double>(last, 0.0);
            curSupport *= localSupport;
          }

          ptrResult[i] += curSupport;
        }
      }
    }

    template<>
    inline void OCLLinear<float>::multDefault(
      sg::base::DataMatrixSP* level,
      sg::base::DataMatrixSP* index,
      sg::base::DataMatrixSP* /*mask*/, //unused for this specialization
      sg::base::DataMatrixSP* /*offset*/, //unused for this specialization
      sg::base::DataMatrixSP* dataset,
      sg::base::DataVectorSP& alpha,
      sg::base::DataVectorSP& result,
      const size_t start_index_grid,
      const size_t end_index_grid,
      const size_t start_index_data,
      const size_t end_index_data) {
      float* ptrLevel = level->getPointer();
      float* ptrIndex = index->getPointer();
      float* ptrAlpha = alpha.getPointer();
      float* ptrData = dataset->getPointer();
      float* ptrResult = result.getPointer();
      size_t result_size = result.getSize();
      size_t dims = dataset->getNrows();

      for (size_t i = start_index_data; i < end_index_data; i++) {
        for (size_t j = start_index_grid; j < end_index_grid; j++) {
          float curSupport = ptrAlpha[j];

          for (size_t d = 0; d < dims; d++) {
            float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
            float index_calc = eval - (ptrIndex[(j * dims) + d]);
            float abs = static_cast<float>(fabs(index_calc));
            float last = 1.0f - abs;
            float localSupport = std::max<float>(last, 0.0f);
            curSupport *= localSupport;
          }

          ptrResult[i] += curSupport;
        }
      }
    }

    template<>
    inline void OCLLinear<double>::multTransposeDefault(
      sg::base::DataMatrix* level,
      sg::base::DataMatrix* index,
      sg::base::DataMatrix* /*mask*/, //unused for this specialization
      sg::base::DataMatrix* /*offset*/, //unused for this specialization
      sg::base::DataMatrix* dataset,
      sg::base::DataVector& source,
      sg::base::DataVector& result,
      const size_t start_index_grid,
      const size_t end_index_grid,
      const size_t start_index_data,
      const size_t end_index_data) {
      double* ptrLevel = level->getPointer();
      double* ptrIndex = index->getPointer();
      double* ptrSource = source.getPointer();
      double* ptrData = dataset->getPointer();
      double* ptrResult = result.getPointer();
      size_t source_size = source.getSize();
      size_t dims = dataset->getNrows();

      for (size_t j = start_index_grid; j < end_index_grid; j++) {
        for (size_t i = start_index_data; i < end_index_data; i++) {
          double curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
            double index_calc = eval - (ptrIndex[(j * dims) + d]);
            double abs = fabs(index_calc);
            double last = 1.0 - abs;
            double localSupport = std::max<double>(last, 0.0);
            curSupport *= localSupport;
          }

          ptrResult[j] += curSupport;
        }
      }
    }

    template<>
    inline void OCLLinear<float>::multTransposeDefault(
      sg::base::DataMatrixSP* level,
      sg::base::DataMatrixSP* index,
      sg::base::DataMatrixSP* /*mask*/, //unused for this specialization
      sg::base::DataMatrixSP* /*offset*/, //unused for this specialization
      sg::base::DataMatrixSP* dataset,
      sg::base::DataVectorSP& source,
      sg::base::DataVectorSP& result,
      const size_t start_index_grid,
      const size_t end_index_grid,
      const size_t start_index_data,
      const size_t end_index_data) {
      float* ptrLevel = level->getPointer();
      float* ptrIndex = index->getPointer();
      float* ptrSource = source.getPointer();
      float* ptrData = dataset->getPointer();
      float* ptrResult = result.getPointer();
      size_t source_size = source.getSize();
      size_t dims = dataset->getNrows();

      for (size_t j = start_index_grid; j < end_index_grid; j++) {
        for (size_t i = start_index_data; i < end_index_data; i++) {
          float curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
            float index_calc = eval - (ptrIndex[(j * dims) + d]);
            float abs = (float)fabs(index_calc);
            float last = 1.0f - abs;
            float localSupport = std::max<float>(last, 0.0f);
            curSupport *= localSupport;
          }

          ptrResult[j] += curSupport;
        }
      }
    }
  }
}
#endif // OCLLINEAR_HPP
