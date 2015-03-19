// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TYPESDATADRIVEN_HPP
#define TYPESDATADRIVEN_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace datadriven {

    /**
     * struct to encapsulate the classsifiers quality by its
     * characteristic numbers
     */
    struct ClassificatorQuality {
      /// number of true positive classified instances
      size_t truePositive_;
      /// number of true negative classified instances
      size_t trueNegative_;
      /// number of false positive classified instances
      size_t falsePositive_;
      /// number of false negative classified instances
      size_t falseNegative_;
    };

    /**
     * strcut to encapsulate the learner's timings
     * during training
     */
    struct LearnerTiming {
      /// complete learning time
      float_t timeComplete_;
      /// time to apply B (including data transfer to eventually used accelerators)
      float_t timeMultComplete_;
      /// pure application time of B
      float_t timeMultCompute_;
      /// time to apply B^T (including data transfer to eventually used accelerators)
      float_t timeMultTransComplete_;
      /// pure application time of B^T
      float_t timeMultTransCompute_;
      /// time of regularization
      float_t timeRegularization_;
      /// number of executed Floating Point operations
      float_t GFlop_;
      /// number of transferred Gbytes
      float_t GByte_;
    };

  }

}

#endif /* TYPESDATADRIVEN_HPP */