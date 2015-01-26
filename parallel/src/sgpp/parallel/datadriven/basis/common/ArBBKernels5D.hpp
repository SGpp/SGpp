/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ARBBKERNELS5D_HPP
#define ARBBKERNELS5D_HPP

#include <cstdlib>

namespace sg {

  namespace parallel {

    class ArBBKernels5D {
      private:
        bool isMultTransSPfirst;
        bool isMultSPfirst;
        bool isMultTransfirst;
        bool isMultfirst;

      public:
        ArBBKernels5D();

        ~ArBBKernels5D();

        double multTransArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

        double multArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims);

        double multTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

        double multSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims);

        void resetKernels();
    };

  }

}

#endif /* ARBBKERNELS5D_HPP */
