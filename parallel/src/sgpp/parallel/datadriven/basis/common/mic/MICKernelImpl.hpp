/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef MICKERNELIMPL_HPP
#define MICKERNELIMPL_HPP

#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

namespace sg {
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

      void uploadGrid(sg::base::DataMatrix* level, sg::base::DataMatrix* index, sg::base::DataMatrix* mask, sg::base::DataMatrix* offset);

      void uploadData(sg::base::DataMatrix* data);

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
