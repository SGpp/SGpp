// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MICKERNELIMPL_HPP
#define MICKERNELIMPL_HPP

#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace mic {
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
      extern double* ptrLevel;
      extern double* ptrIndex;
      extern double* ptrMask;
      extern double* ptrOffset;
      extern double* ptrData;
      extern double* ptrDataMic;
      extern double* ptrAlphaMic;
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
      extern int number_mic_devices;
      extern bool multicard_multtrans_fast;
      extern double** tempgrid;

      void uploadGrid(SGPP::base::DataMatrix* level, SGPP::base::DataMatrix* index, SGPP::base::DataMatrix* mask, SGPP::base::DataMatrix* offset);

      void uploadData(SGPP::base::DataMatrix* data);

      void deleteGrid();

      void deleteData();

      /**
       * @brief transferResultMult workaround for intel compiler bug: it's not possible to have an into clause in a templated class
       * @param offset offset of the result to copy
       * @param size size of the result to copy
       * @param device which device do we want to copy the result from
       * @param ptrResult array on the host for the result
       */
      void transferResultMult(size_t offset, size_t size, size_t device, double* ptrResult);

      /**
       * @brief transferResultMultTrans workaround for intel compiler bug: it's not possible to have an into clause in a templated class
       * @param offset offset of the result to copy
       * @param size size of the result to copy
       * @param device which device do we want to copy the result from
       * @param ptrResult array on the host for the result
       */
      void transferResultMultTrans(size_t offset, size_t size, size_t device, double* ptrResult);

      void transferInputMult(size_t offsetAlpha, size_t chunkAlpha, double* ptrAlpha,
                             size_t offsetResult, size_t chunkResult, double* ptrResult, size_t device);
      void transferInputMultTrans(size_t offsetSource, size_t chunkSource, double* ptrSource,
                                  size_t offsetResult, size_t chunkResult, double* ptrResult, size_t device);
    }
  }
}
#endif // MICKERNELIMPL_HPP