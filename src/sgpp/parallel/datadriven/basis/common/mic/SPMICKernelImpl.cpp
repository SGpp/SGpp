/* ****************************************************************************
* Copyright (C) 2010-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifdef USEMIC

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <iostream>
#include <sstream>
#include <cstring>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include "parallel/datadriven/basis/common/mic/SPMICKernelImpl.hpp"
#include <offload.h>


namespace sg {
  namespace parallel {
    namespace micsp {
      float* ptrLevel = NULL;
      float* ptrIndex = NULL;
      float* ptrMask = NULL;
      float* ptrOffset = NULL;
      float* ptrData = NULL;
      float* ptrDataMic = NULL;
      float* ptrAlphaMic = NULL;

      int init_number_mic_devices() {
#ifdef __INTEL_OFFLOAD
        const char* num_mic_devices_env = getenv("SGPP_NUM_MIC_DEVICES");
        int num_mic_devs = _Offload_number_of_devices();

        if (num_mic_devices_env != NULL) {
          int num_mic_devices_limit = (int)(strtoul (num_mic_devices_env, NULL, 0));

          if (num_mic_devices_limit != 0) {
            num_mic_devs = std::min<int>(num_mic_devs, num_mic_devices_limit);
          } else {
            std::cout << "Ignoring value: \"" << num_mic_devices_env << "\" for SGPP_NUM_MIC_DEVICES" << std::endl;
          }
        }

        if (num_mic_devs == -1) {
          std::cout << "Warning: there are no MIC devices available on this machine, execution will be slow." << std::endl;
        } else {
          std::cout << "Using " << num_mic_devs << " MIC devices ()." << std::endl;
        }

        return num_mic_devs;
#else
        return 0;
#endif
      }

      int number_mic_devices = init_number_mic_devices();
      extern float** tempgrid = NULL;

      bool init_multicard_multtrans_fast() {
        const char* multicard_multtrans_partitioning = getenv("SGPP_MULTICARD_MIC"); // =fast|safe

        if (multicard_multtrans_partitioning == NULL) {
          multicard_multtrans_partitioning = "fast";
        }

        std::cout << "SGPP_MULTICARD_MIC: " << multicard_multtrans_partitioning << std::endl;

        if (strcmp(multicard_multtrans_partitioning, "safe") == 0) {
          return false;
        } else if (strcmp(multicard_multtrans_partitioning, "fast") == 0) {
          return true;
        } else {
          std::cerr << "SGPP_MULTICARD_MIC must be either 'fast' or 'safe'. Assuming 'safe'" << std::endl;
          return true;
        }
      }

      bool multicard_multtrans_fast = init_multicard_multtrans_fast();

      void uploadGrid(sg::base::DataMatrixSP* level, sg::base::DataMatrixSP* index, sg::base::DataMatrixSP* mask, sg::base::DataMatrixSP* offset) {
        size_t storageSize = level->getNrows();
        size_t dims = level->getNcols();

        if (level != NULL) {
          ptrLevel = level->getPointer();
#ifdef __INTEL_OFFLOAD
          #pragma omp parallel for schedule(static,1)

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) in(ptrLevel:length(storageSize*dims) free_if(0) alloc_if(1) align(64))
          }

#endif
        }

        if (index != NULL) {
          ptrIndex = index->getPointer();
#ifdef __INTEL_OFFLOAD
          #pragma omp parallel for schedule(static,1)

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) in(ptrIndex:length(storageSize*dims) free_if(0) alloc_if(1) align(64))
          }

#endif
        }

        if (mask != NULL) {
          ptrMask = mask->getPointer();
#ifdef __INTEL_OFFLOAD
          #pragma omp parallel for schedule(static,1)

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) in(ptrMask:length(storageSize*dims) free_if(0) alloc_if(1) align(64))
          }

#endif
        }

        if (offset != NULL) {
          ptrOffset = offset->getPointer();
#ifdef __INTEL_OFFLOAD
          #pragma omp parallel for schedule(static,1)

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) in(ptrOffset:length(storageSize*dims) free_if(0) alloc_if(1) align(64))
          }

#endif
        }

#ifdef __INTEL_OFFLOAD
        ptrAlphaMic =  new float[storageSize];
        memset(ptrAlphaMic, 0, sizeof(float)*storageSize);
        #pragma omp parallel for schedule(static,1)

        for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) in(ptrAlphaMic:length(storageSize) free_if(0) alloc_if(1) align(64))
        }

#endif
      }

      void uploadData(sg::base::DataMatrixSP* data) {
        ptrData = data->getPointer();
        size_t datasize = data->getNcols();
        size_t dims = data->getNrows();
#ifdef __INTEL_OFFLOAD
        ptrDataMic = new float[datasize];
        memset(ptrDataMic, 0, sizeof(float)*datasize);

        for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) \
  in(ptrData:length(datasize*dims) free_if(0) alloc_if(1) align(64)) \
  in(ptrDataMic:length(datasize) free_if(0) alloc_if(1) align(64))
        }

#endif
      }

      void deleteGrid() {
        if (ptrLevel != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrLevel:length(0) alloc_if(0) free_if(1)) // alloc_if(0) is default for nocopy, length is ignored for free
          }

#endif

          ptrLevel = NULL;
        }

        if (ptrIndex != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrIndex:length(0) alloc_if(0) free_if(1)) // alloc_if(0) is default for nocopy, length is ignored for free
          }

#endif

          ptrIndex = NULL;
        }

        if (ptrMask != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrMask:length(0) alloc_if(0) free_if(1)) // alloc_if(0) is default for nocopy, length is ignored for free
          }

#endif

          ptrMask = NULL;
        }

        if (ptrOffset != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrOffset:length(0) alloc_if(0) free_if(1)) // alloc_if(0) is default for nocopy, length is ignored for free
          }

#endif

          ptrOffset = NULL;
        }

        if (ptrAlphaMic != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrAlphaMic:length(0) alloc_if(0) free_if(1)) // alloc_if(0) is default for nocopy, length is ignored for free
          }

#endif

          delete[] ptrAlphaMic;
          ptrAlphaMic = NULL;
        }
      }

      void deleteData() {
        if (ptrData != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrData:length(0) alloc_if(0) free_if(1))
          }

#endif

          ptrData = NULL;
        }

        if (ptrDataMic != NULL) {
#ifdef __INTEL_OFFLOAD

          for (size_t d = 0; d < number_mic_devices; d++) {
#pragma offload_transfer target(mic:d) nocopy(ptrDataMic:length(0) alloc_if(0) free_if(1))
          }

#endif

          delete[] ptrDataMic;
          ptrDataMic = NULL;
        }
      }

#ifdef __INTEL_OFFLOAD
      /**
       * @brief transferResultMult workaround for intel compiler bug: it's not possible to have an into clause in a templated class
       * @param offset offset of the result to copy
       * @param size size of the result to copy
       * @param device which device do we want to copy the result from
       * @param ptrResult array on the host for the result
       */
      void transferResultMult(size_t offset, size_t size, size_t device, float* ptrResult) {
#pragma offload_transfer target(mic:device) out(ptrDataMic[offset:size] : into(ptrResult[offset:size]) free_if(0) alloc_if(0) align(64))
      }

      /**
       * @brief transferResultMultTrans workaround for intel compiler bug: it's not possible to have an into clause in a templated class
       * @param offset offset of the result to copy
       * @param size size of the result to copy
       * @param device which device do we want to copy the result from
       * @param ptrResult array on the host for the result
       */
      void transferResultMultTrans(size_t offset, size_t size, size_t device, float* ptrResult) {
#pragma offload_transfer target(mic:device) out(ptrAlphaMic[offset:size] : into(ptrResult[offset:size]) free_if(0) alloc_if(0) align(64))
      }

      void transferInputMult(size_t offsetAlpha, size_t chunkAlpha, float* ptrAlpha,
                             size_t offsetResult, size_t chunkResult, float* ptrResult, size_t device) {
#pragma offload_transfer target(mic:device) in(ptrAlpha[offsetAlpha:chunkAlpha] : into(ptrAlphaMic[offsetAlpha:chunkAlpha]) free_if(0) alloc_if(0) align(64))
#pragma offload_transfer target(mic:device) in(ptrResult[offsetResult:chunkResult] : into(ptrDataMic[offsetResult:chunkResult]) free_if(0) alloc_if(0) align(64))
      }
      void transferInputMultTrans(size_t offsetSource, size_t chunkSource, float* ptrSource,
                                  size_t offsetResult, size_t chunkResult, float* ptrResult, size_t device) {
#pragma offload_transfer target(mic:device) in(ptrSource[offsetSource:chunkSource] : into(ptrDataMic[offsetSource:chunkSource]) free_if(0) alloc_if(0) align(64))
#pragma offload_transfer target(mic:device) in(ptrResult[offsetResult:chunkResult] : into(ptrAlphaMic[offsetResult:chunkResult]) free_if(0) alloc_if(0) align(64))
      }
#endif

    }
  }
}

#endif // USEMIC
