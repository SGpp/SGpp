// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ARBBKERNELS4D_HPP
#define ARBBKERNELS4D_HPP

#include <cstdlib>

#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace parallel {

class ArBBKernels4D {
 private:
  bool isMultTransSPfirst;
  bool isMultSPfirst;
  bool isMultTransfirst;
  bool isMultfirst;

 public:
  ArBBKernels4D();

  ~ArBBKernels4D();

  double multTransArBB(double* ptrSource, double* ptrData, double* ptrLevel,
                       double* ptrIndex, double* ptrGlobalResult, size_t sourceSize,
                       size_t storageSize, size_t dims);

  double multArBB(double* ptrAlpha, double* ptrData, double* ptrLevel,
                  double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize,
                  size_t dims);

  double multTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel,
                         float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize,
                         size_t dims);

  double multSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel,
                    float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize,
                    size_t dims);

  void resetKernels();
};

}

}

#endif /* ARBBKERNELS4D_HPP */