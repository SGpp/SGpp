/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"

namespace sg {
  namespace parallel {

    OperationMultipleEvalVectorizedSP::OperationMultipleEvalVectorizedSP(base::GridStorage* storage, base::DataMatrixSP* dataset) {
      this->storage_ = storage;
      this->dataset_ = dataset;
      this->level_ = NULL;
      this->index_ = NULL;
      this->mask_ = NULL;
      this->offset_ = NULL;
      this->myTimer_ = new sg::base::SGppStopwatch();
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
