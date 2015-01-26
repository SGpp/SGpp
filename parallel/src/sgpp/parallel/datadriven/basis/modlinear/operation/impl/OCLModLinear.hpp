/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OCLMODLINEAR_HPP
#define OCLMODLINEAR_HPP

#include <sstream>
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernelBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    template<typename real_type>
    class OCLModLinear : public OCLKernelBase {
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
          stream_program_src << "	int globalIdx = get_global_id(0);" << std::endl;
          stream_program_src << "	int localIdx = get_local_id(0);" << std::endl;
          stream_program_src << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << "	__local " << getType<real_type>::asString() << " locLevel[" << dims* local_workgroup_size << "];" << std::endl;
          stream_program_src << "	__local " << getType<real_type>::asString() << " locIndex[" << dims* local_workgroup_size << "];" << std::endl;
          stream_program_src << "	__local " << getType<real_type>::asString() << " locAlpha[" << local_workgroup_size << "];" << std::endl;
          stream_program_src << std::endl;
#endif
          stream_program_src << "	" << getType<real_type>::asString() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
          stream_program_src << "	" << getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl << std::endl;
          stream_program_src << "	// Create registers for the data" << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "	" << getType<real_type>::asString() << " data_" << d << " = ptrData[globalIdx+(resultSize*" << d << ")];" << std::endl;
          }

          stream_program_src << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << "	// Iterate over all grid points (fast ones, with cache)" << std::endl;
          stream_program_src << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
          stream_program_src << " uint fastChunkSizeGrid = (chunkSizeGrid / " << local_workgroup_size << ") * " << local_workgroup_size << ";" << std::endl;
          stream_program_src << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+=" << local_workgroup_size << ")" << std::endl;
          stream_program_src << "	{" << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "		locLevel[(localIdx*" << dims << ")+" << d << "] = ptrLevel[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
            stream_program_src << "		locIndex[(localIdx*" << dims << ")+" << d << "] = ptrIndex[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
          }

          stream_program_src << "		locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
          stream_program_src << "		barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          stream_program_src << std::endl;
          stream_program_src << "		for(int k = 0; k < " << local_workgroup_size << "; k++)" << std::endl;
          stream_program_src << "		{" << std::endl;
          stream_program_src << "			curSupport = locAlpha[k];" << std::endl << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "			if ((locLevel[(k*" << dims << ")+" << d << "]) == 2.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "			{" << std::endl;
            stream_program_src << "				curSupport *= 1.0" << getType<real_type>::constSuffix() << ";" << std::endl;
            stream_program_src << "			}" << std::endl;
            stream_program_src << "			else if ((locIndex[(k*" << dims << ")+" << d << "]) == 1.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "			{" << std::endl;
            stream_program_src << "				curSupport *= max(2.0" << getType<real_type>::constSuffix() << " - ( (locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << ") ), 0.0" << getType<real_type>::constSuffix() << ") ;" << std::endl;
            stream_program_src << "			}" << std::endl;
            stream_program_src << "			else if ((locIndex[(k*" << dims << ")+" << d << "]) == ((locLevel[(k*" << dims << ")+" << d << "]) - 1.0" << getType<real_type>::constSuffix() << ") )" << std::endl;
            stream_program_src << "			{" << std::endl;
            stream_program_src << "				curSupport *= max(( (locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (locIndex[(k*" << dims << ")+" << d << "]) + 1.0" << getType<real_type>::constSuffix() << ", 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "			}" << std::endl;
            stream_program_src << "			else " << std::endl;
            stream_program_src << "			{" << std::endl;
            stream_program_src << "				curSupport *= max(1.0" << getType<real_type>::constSuffix() << " - fabs( ( (locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (locIndex[(k*" << dims << ")+" << d << "]) ), 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "			}" << std::endl;
          }

          stream_program_src << "			myResult += curSupport;" << std::endl;
          stream_program_src << "		}" << std::endl;
          stream_program_src << std::endl;
          stream_program_src << "		barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          stream_program_src << "	}" << std::endl;
          stream_program_src << std::endl;
          stream_program_src << "	// Iterate over all grid points (slow ones, without cache)" << std::endl;
          stream_program_src << " for(int m = start_grid + fastChunkSizeGrid; m < end_grid; m++)" << std::endl;
          stream_program_src << "	{" << std::endl;
          stream_program_src << "		curSupport = ptrAlpha[m];" << std::endl << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "		if ((ptrLevel[(m*" << dims << ")+" << d << "]) == 2.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= 1.0" << getType<real_type>::constSuffix() << ";" << std::endl;
            stream_program_src << "		}" << std::endl;
            stream_program_src << "		else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == 1.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= max(2.0" << getType<real_type>::constSuffix() << " - ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ), 0.0" << getType<real_type>::constSuffix() << ") ;" << std::endl;
            stream_program_src << "		}" << std::endl;
            stream_program_src << "		else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == ((ptrLevel[(m*" << dims << ")+" << d << "]) - 1.0" << getType<real_type>::constSuffix() << ") )" << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= max(( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (ptrIndex[(m*" << dims << ")+" << d << "]) + 1.0" << getType<real_type>::constSuffix() << ", 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "		}" << std::endl;
            stream_program_src << "		else " << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= max(1.0" << getType<real_type>::constSuffix() << " - fabs( ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (ptrIndex[(m*" << dims << ")+" << d << "]) ), 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "		}" << std::endl;
          }

          stream_program_src << "		myResult += curSupport;" << std::endl;
          stream_program_src << "	}" << std::endl;
#else
          stream_program_src << "	// Iterate over all grid points (slow ones, without cache)" << std::endl;
          stream_program_src << " for(int m = start_grid; m < end_grid; m++)" << std::endl;
          stream_program_src << "	{" << std::endl;
          stream_program_src << "		curSupport = ptrAlpha[m];" << std::endl << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "		if ((ptrLevel[(m*" << dims << ")+" << d << "]) == 2.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= 1.0" << getType<real_type>::constSuffix() << ";" << std::endl;
            stream_program_src << "		}" << std::endl;
            stream_program_src << "		else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == 1.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= max(2.0" << getType<real_type>::constSuffix() << " - ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ), 0.0" << getType<real_type>::constSuffix() << ") ;" << std::endl;
            stream_program_src << "		}" << std::endl;
            stream_program_src << "		else if ((ptrIndex[(m*" << dims << ")+" << d << "]) == ((ptrLevel[(m*" << dims << ")+" << d << "]) - 1.0" << getType<real_type>::constSuffix() << ") )" << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= max(( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (ptrIndex[(m*" << dims << ")+" << d << "]) + 1.0" << getType<real_type>::constSuffix() << ", 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "		}" << std::endl;
            stream_program_src << "		else " << std::endl;
            stream_program_src << "		{" << std::endl;
            stream_program_src << "			curSupport *= max(1.0" << getType<real_type>::constSuffix() << " - fabs( ( (ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << ") ) - (ptrIndex[(m*" << dims << ")+" << d << "]) ), 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
            stream_program_src << "		}" << std::endl;
          }

          stream_program_src << "		myResult += curSupport;" << std::endl;
          stream_program_src << "	}" << std::endl;
#endif
          stream_program_src << std::endl;
          stream_program_src << "	ptrResult[globalIdx] = myResult;" << std::endl;
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
          stream_program_src << "	int globalIdx = get_global_id(0);" << std::endl;
          stream_program_src << "	int localIdx = get_local_id(0);" << std::endl;
          stream_program_src << std::endl;
          stream_program_src << "	" << getType<real_type>::asString() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
          stream_program_src << "	" << getType<real_type>::asString() << " myResult = ptrResult[globalIdx];" << std::endl << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << "	__local " << getType<real_type>::asString() << " locData[" << dims* local_workgroup_size << "];" << std::endl;
          stream_program_src << "	__local " << getType<real_type>::asString() << " locSource[" << local_workgroup_size << "];" << std::endl << std::endl;
#endif

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "	" << getType<real_type>::asString() << " level_" << d << " = ptrLevel[(globalIdx*" << dims << ")+" << d << "];" << std::endl;
            stream_program_src << "	" << getType<real_type>::asString() << " index_" << d << " = ptrIndex[(globalIdx*" << dims << ")+" << d << "];" << std::endl;
          }

          stream_program_src << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << "	// Iterate over all grid points" << std::endl;
          stream_program_src << " for(int i = start_data; i < end_data; i+=" << local_workgroup_size << ")" << std::endl;
          stream_program_src << "	{" << std::endl;

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "		locData[(" << d << "*" << local_workgroup_size << ")+(localIdx)] = ptrData[(" << d << "*sourceSize)+(localIdx+i)];" << std::endl;
          }

          stream_program_src << "		locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
          stream_program_src << "		barrier(CLK_LOCAL_MEM_FENCE);" << std::endl << std::endl;
          stream_program_src << "		for(int k = 0; k < " << local_workgroup_size << "; k++)" << std::endl;
          stream_program_src << "		{" << std::endl;
          stream_program_src << "			curSupport = locSource[k];" << std::endl << std::endl;
#else
          stream_program_src << "   for(int k = start_data; k < end_data; k++)" << std::endl;
          stream_program_src << "		{" << std::endl;
          stream_program_src << "			curSupport = ptrSource[k];" << std::endl << std::endl;
#endif

          for (size_t d = 0; d < dims; d++) {
            stream_program_src << "			if ((level_" << d << ") == 2.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "			{" << std::endl;
            stream_program_src << "				curSupport *= 1.0" << getType<real_type>::constSuffix() << ";" << std::endl;
            stream_program_src << "			}" << std::endl;
            stream_program_src << "			else if ((index_" << d << ") == 1.0" << getType<real_type>::constSuffix() << ")" << std::endl;
            stream_program_src << "			{" << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
            stream_program_src << "				curSupport *= max(2.0" << getType<real_type>::constSuffix() << " - ( (level_" << d << ") * (locData[(" << d << "*" << local_workgroup_size << ")+k]) ), 0.0" << getType<real_type>::constSuffix() << ") ;" << std::endl;
#else
            stream_program_src << "				curSupport *= max(2.0" << getType<real_type>::constSuffix() << " - ( (level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]) ), 0.0" << getType<real_type>::constSuffix() << ") ;" << std::endl;
#endif
            stream_program_src << "			}" << std::endl;
            stream_program_src << "			else if ((index_" << d << ") == ((level_" << d << ") - 1.0" << getType<real_type>::constSuffix() << ") )" << std::endl;
            stream_program_src << "			{" << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
            stream_program_src << "				curSupport *= max(( (level_" << d << ") * (locData[(" << d << "*" << local_workgroup_size << ")+k]) ) - (index_" << d << ") + 1.0" << getType<real_type>::constSuffix() << ", 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
#else
            stream_program_src << "				curSupport *= max(( (level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]) ) - (index_" << d << ") + 1.0, 0.0);" << std::endl;
#endif
            stream_program_src << "			}" << std::endl;
            stream_program_src << "			else " << std::endl;
            stream_program_src << "			{" << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
            stream_program_src << "				curSupport *= max(1.0" << getType<real_type>::constSuffix() << " - fabs( ( (level_" << d << ") * (locData[(" << d << "*" << local_workgroup_size << ")+k]) ) - (index_" << d << ") ), 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
#else
            stream_program_src << "				curSupport *= max(1.0" << getType<real_type>::constSuffix() << " - fabs( ( (level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]) ) - (index_" << d << ") ), 0.0" << getType<real_type>::constSuffix() << ");" << std::endl;
#endif
            stream_program_src << "			}" << std::endl;
          }

          stream_program_src << std::endl << "		myResult += curSupport;" << std::endl;
          stream_program_src << "		}" << std::endl << std::endl;
#ifdef USEOCL_LOCAL_MEMORY
          stream_program_src << "		barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          stream_program_src << "	}" << std::endl;
#endif
          stream_program_src << "	ptrResult[globalIdx] = myResult;" << std::endl;
          stream_program_src << "}" << std::endl;

          return stream_program_src.str();
        }
      public:
        static inline void multDefault(
          SGPP::base::DataMatrix* level,
          SGPP::base::DataMatrix* index,
          SGPP::base::DataMatrix* /*mask*/, //unused for this specialization
          SGPP::base::DataMatrix* /*offset*/, //unused for this specialization
          SGPP::base::DataMatrix* dataset,
          SGPP::base::DataVector& alpha,
          SGPP::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
        static inline void multDefault(
          SGPP::base::DataMatrixSP* level,
          SGPP::base::DataMatrixSP* index,
          SGPP::base::DataMatrixSP* /*mask*/, //unused for this specialization
          SGPP::base::DataMatrixSP* /*offset*/, //unused for this specialization
          SGPP::base::DataMatrixSP* dataset,
          SGPP::base::DataVectorSP& alpha,
          SGPP::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
        static inline void multTransposeDefault(
          SGPP::base::DataMatrix* level,
          SGPP::base::DataMatrix* index,
          SGPP::base::DataMatrix* /*mask*/, //unused for this specialization
          SGPP::base::DataMatrix* /*offset*/, //unused for this specialization
          SGPP::base::DataMatrix* dataset,
          SGPP::base::DataVector& source,
          SGPP::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
        static inline void multTransposeDefault(
          SGPP::base::DataMatrixSP* level,
          SGPP::base::DataMatrixSP* index,
          SGPP::base::DataMatrixSP* /*mask*/, //unused for this specialization
          SGPP::base::DataMatrixSP* /*offset*/, //unused for this specialization
          SGPP::base::DataMatrixSP* dataset,
          SGPP::base::DataVectorSP& source,
          SGPP::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data);
    };
    template<>
    inline void OCLModLinear<double>::multDefault(
      SGPP::base::DataMatrix* level,
      SGPP::base::DataMatrix* index,
      SGPP::base::DataMatrix* /*mask*/, //unused for this specialization
      SGPP::base::DataMatrix* /*offset*/, //unused for this specialization
      SGPP::base::DataMatrix* dataset,
      SGPP::base::DataVector& alpha,
      SGPP::base::DataVector& result,
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
            if (ptrLevel[(j * dims) + d] == 2.0) {
              // nothing to do (mult with 1)
            } else if (ptrIndex[(j * dims) + d] == 1.0) {
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              eval = 2.0 - eval;
              double localSupport = std::max<double>(eval, 0.0);
              curSupport *= localSupport;
            } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0)) {
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              double index_calc = eval - (ptrIndex[(j * dims) + d]);
              double last = 1.0 + index_calc;
              double localSupport = std::max<double>(last, 0.0);
              curSupport *= localSupport;
            } else {
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              double index_calc = eval - (ptrIndex[(j * dims) + d]);
              double abs = fabs(index_calc);
              double last = 1.0 - abs;
              double localSupport = std::max<double>(last, 0.0);
              curSupport *= localSupport;
            }
          }

          ptrResult[i] += curSupport;
        }
      }
    }

    template<>
    inline void OCLModLinear<float>::multDefault(
      SGPP::base::DataMatrixSP* level,
      SGPP::base::DataMatrixSP* index,
      SGPP::base::DataMatrixSP* /*mask*/, //unused for this specialization
      SGPP::base::DataMatrixSP* /*offset*/, //unused for this specialization
      SGPP::base::DataMatrixSP* dataset,
      SGPP::base::DataVectorSP& alpha,
      SGPP::base::DataVectorSP& result,
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
            if (ptrLevel[(j * dims) + d] == 2.0f) {
              // nothing to do (mult with 1)
            } else if (ptrIndex[(j * dims) + d] == 1.0f) {
              float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              eval = 2.0f - eval;
              float localSupport = std::max<float>(eval, 0.0f);
              curSupport *= localSupport;
            } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
              float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              float index_calc = eval - (ptrIndex[(j * dims) + d]);
              float last = 1.0f + index_calc;
              float localSupport = std::max<float>(last, 0.0f);
              curSupport *= localSupport;
            } else {
              float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              float index_calc = eval - (ptrIndex[(j * dims) + d]);
              float abs = static_cast<float>(fabs(index_calc));
              float last = 1.0f - abs;
              float localSupport = std::max<float>(last, 0.0f);
              curSupport *= localSupport;
            }
          }

          ptrResult[i] += curSupport;
        }
      }
    }


    template<>
    inline void OCLModLinear<double>::multTransposeDefault(
      SGPP::base::DataMatrix* level,
      SGPP::base::DataMatrix* index,
      SGPP::base::DataMatrix* /*mask*/, //unused for this specialization
      SGPP::base::DataMatrix* /*offset*/, //unused for this specialization
      SGPP::base::DataMatrix* dataset,
      SGPP::base::DataVector& source,
      SGPP::base::DataVector& result,
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
            if (ptrLevel[(j * dims) + d] == 2.0) {
              curSupport *= 1.0;
            } else if (ptrIndex[(j * dims) + d] == 1.0) {
              curSupport *= std::max<double>(2.0 - ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])), 0.0);
            } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0)) {
              curSupport *= std::max<double>(((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])) - ptrIndex[(j * dims) + d] + 1.0, 0.0);
            } else {
              curSupport *= std::max<double>(1.0 - fabs( ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])) - ptrIndex[(j * dims) + d] ), 0.0);
            }
          }

          ptrResult[j] += curSupport;
        }
      }

    }

    template<>
    inline void OCLModLinear<float>::multTransposeDefault(
      SGPP::base::DataMatrixSP* level,
      SGPP::base::DataMatrixSP* index,
      SGPP::base::DataMatrixSP* /*mask*/, //unused for this specialization
      SGPP::base::DataMatrixSP* /*offset*/, //unused for this specialization
      SGPP::base::DataMatrixSP* dataset,
      SGPP::base::DataVectorSP& source,
      SGPP::base::DataVectorSP& result,
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
            if (ptrLevel[(j * dims) + d] == 2.0f) {
              curSupport *= 1.0f;
            } else if (ptrIndex[(j * dims) + d] == 1.0f) {
              curSupport *= std::max<float>(2.0f - ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])), 0.0f);
            } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
              curSupport *= std::max<float>(((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])) - ptrIndex[(j * dims) + d] + 1.0f, 0.0f);
            } else {
              curSupport *= std::max<float>(1.0f - (float)fabs( ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])) - ptrIndex[(j * dims) + d] ), 0.0f);
            }
          }

          ptrResult[j] += curSupport;
        }
      }
    }
  }
}
#endif // OCLMODLINEAR_HPP
