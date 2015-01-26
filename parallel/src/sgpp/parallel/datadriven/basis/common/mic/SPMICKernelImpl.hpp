/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPMICKERNELIMPL_HPP
#define SPMICKERNELIMPL_HPP

#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

namespace sg {
  namespace parallel {
    namespace micsp {
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
      extern float* ptrLevel;
      extern float* ptrIndex;
      extern float* ptrMask;
      extern float* ptrOffset;
      extern float* ptrData;
      extern float* ptrDataMic;
      extern float* ptrAlphaMic;
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
      extern int number_mic_devices;
      extern bool multicard_multtrans_fast;
      extern float** tempgrid;

      void uploadGrid(sg::base::DataMatrixSP* level, sg::base::DataMatrixSP* index, sg::base::DataMatrixSP* mask, sg::base::DataMatrixSP* offset);

      void uploadData(sg::base::DataMatrixSP* data);

      void deleteGrid();

      void deleteData();

      /**
       * @brief transferResultMult workaround for intel compiler bug: it's not possible to have an into clause in a templated class
       * @param offset offset of the result to copy
       * @param size size of the result to copy
       * @param device which device do we want to copy the result from
       * @param ptrResult array on the host for the result
       */
      void transferResultMult(size_t offset, size_t size, size_t device, float* ptrResult);

      /**
       * @brief transferResultMultTrans workaround for intel compiler bug: it's not possible to have an into clause in a templated class
       * @param offset offset of the result to copy
       * @param size size of the result to copy
       * @param device which device do we want to copy the result from
       * @param ptrResult array on the host for the result
       */
      void transferResultMultTrans(size_t offset, size_t size, size_t device, float* ptrResult);

      void transferInputMult(size_t offsetAlpha, size_t chunkAlpha, float* ptrAlpha,
                             size_t offsetResult, size_t chunkResult, float* ptrResult, size_t device);
      void transferInputMultTrans(size_t offsetSource, size_t chunkSource, float* ptrSource,
                                  size_t offsetResult, size_t chunkResult, float* ptrResult, size_t device);
    }
  }
}
#endif // SPMICKERNELIMPL_HPP
