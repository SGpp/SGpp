// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    enum class OperationMultipleEvalType {
      DEFAULT, STREAMING, SUBSPACELINEAR, ADAPTIVE
    };

    enum class OperationMultipleEvalSubType {
      DEFAULT, SIMPLE, COMBINED, OCL, OCLFAST, OCLFASTMULTIPLATFORM
    };

    class OperationMultipleEvalConfiguration {
      public:
        OperationMultipleEvalType type = OperationMultipleEvalType::DEFAULT;
        OperationMultipleEvalSubType subType = OperationMultipleEvalSubType::DEFAULT;

        //operational - can be set for easier reporting
        std::string name = "DEFAULT";
    };

  }
}
