/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TYPESPARALLEL_HPP
#define TYPESPARALLEL_HPP

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace parallel {

    enum VectorizationType {
      X86SIMD,
      OpenCL,
      Hybrid_X86SIMD_OpenCL,
      MIC,
      Hybrid_X86SIMD_MIC,
      CUDA,
      ArBB
    };

    enum MPIType {
      MPINone,
      MPIAllreduce,
      MPIAsync,
      MPIOnesided,
      MPIAlltoallv,
      MPITrueAsync,
      MPIBigdata
    };

  }

}

#endif /* TYPESPARALLEL_HPP */
