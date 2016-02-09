// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/datadriven/basis/common/mic/MICKernelBase.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernelImplBase.hpp>
#include <sgpp/parallel/datadriven/basis/common/X86SimdKernelBase.hpp>
#include <sgpp/parallel/datadriven/basis/common/SPX86SimdKernelBase.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace parallel {

size_t DMVectorizationPaddingAssistant::getVecWidth(VectorizationType&
    vecType) {
  if (vecType == X86SIMD) {
#ifdef X86_MIC_SYMMETRIC
    VectorizationType t = MIC;
    return getVecWidth(t);
#else
    return X86SimdKernelBase::getChunkDataPoints();
#endif
  }

#ifdef USEOCL
  else if (vecType == OpenCL) {
    return OCLKernelImplBase::getChunkDataPoints();
  } else if (vecType == Hybrid_X86SIMD_OpenCL) {
    return OCLKernelImplBase::getChunkDataPoints();
  }

#endif
  else if (vecType == ArBB) {
    return 16;
  }

#ifdef USEMIC
  else if (vecType == MIC) {
    return MICKernelBase::getChunkDataPoints();
  } else if (vecType == Hybrid_X86SIMD_MIC) {
    return MICKernelBase::getChunkDataPoints();
  }

#endif
  else {
    throw SGPP::base::operation_exception("DMVectorizationPaddingAssistant::getVecWidth : un-supported vector extension!");
  }

  return 0;
}

size_t DMVectorizationPaddingAssistant::getVecWidthSP(VectorizationType&
    vecType) {
  if (vecType == X86SIMD) {
#ifdef X86_MIC_SYMMETRIC
    VectorizationType t = MIC;
    return getVecWidthSP(t);
#else
    return SPX86SimdKernelBase::getChunkDataPoints();
#endif
  }

#ifdef USEOCL
  else if (vecType == OpenCL) {
    return OCLKernelImplBase::getChunkDataPoints();
  } else if (vecType == Hybrid_X86SIMD_OpenCL) {
    return OCLKernelImplBase::getChunkDataPoints();
  }

#endif
  else if (vecType == ArBB) {
    return 16;
  }

#ifdef USECUDA
  else if (vecType == CUDA) {
    return 64;
  }

#endif

#ifdef USEMIC
  else if (vecType == MIC) {
    return SPMICKernelBase::getChunkDataPoints();
  } else if (vecType == Hybrid_X86SIMD_MIC) {
    return SPMICKernelBase::getChunkDataPoints();
  }

#endif
  else {
    throw SGPP::base::operation_exception("DMVectorizationPaddingAssistant::getVecWidthSP : un-supported vector extension!");
  }

  return 0;
}

size_t DMVectorizationPaddingAssistant::padDataset(SGPP::base::DataMatrix&
    dataset, VectorizationType& vecType) {
  size_t vecWidth = getVecWidth(vecType);

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = dataset.getNrows() % vecWidth;
  size_t loopCount = vecWidth - remainder;

  if (loopCount != vecWidth) {
    SGPP::base::DataVector lastRow(dataset.getNcols());
    size_t oldSize = dataset.getNrows();
    dataset.getRow(dataset.getNrows() - 1, lastRow);
    dataset.resize(dataset.getNrows() + loopCount);

    for (size_t i = 0; i < loopCount; i++) {
      dataset.setRow(oldSize + i, lastRow);
    }
  }

  return dataset.getNrows();
}

size_t DMVectorizationPaddingAssistant::padDataset(SGPP::base::DataMatrixSP&
    dataset, VectorizationType vecType) {
  size_t vecWidth = getVecWidthSP(vecType);

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = dataset.getNrows() % vecWidth;
  size_t loopCount = vecWidth - remainder;

  if (loopCount != vecWidth) {
    SGPP::base::DataVectorSP lastRow(dataset.getNcols());
    size_t oldSize = dataset.getNrows();
    dataset.getRow(dataset.getNrows() - 1, lastRow);
    dataset.resize(dataset.getNrows() + loopCount);

    for (size_t i = 0; i < loopCount; i++) {
      dataset.setRow(oldSize + i, lastRow);
    }
  }

  return dataset.getNrows();
}
}

}