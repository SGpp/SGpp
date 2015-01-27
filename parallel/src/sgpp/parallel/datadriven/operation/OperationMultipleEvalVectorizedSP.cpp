// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp>

#include <sgpp/globaldef.hpp>

// TODO David
#if USE_DOUBLE_PRECISION==0

namespace SGPP {
  namespace parallel {

    OperationMultipleEvalVectorizedSP::OperationMultipleEvalVectorizedSP(base::GridStorage* storage, base::DataMatrixSP* dataset) {
      this->storage_ = storage;
      this->dataset_ = dataset;
      this->level_ = NULL;
      this->index_ = NULL;
      this->mask_ = NULL;
      this->offset_ = NULL;
      this->myTimer_ = new SGPP::base::SGppStopwatch();
    }

    OperationMultipleEvalVectorizedSP::~OperationMultipleEvalVectorizedSP() {
      delete myTimer_;

      if (this->level_ != NULL)
        delete this->level_;

      if (this->index_ != NULL)
        delete this->index_;

      if (this->mask_ != NULL)
        delete this->mask_;

      if (this->offset_ != NULL)
        delete this->offset_;
    }

  }
}

#endif
