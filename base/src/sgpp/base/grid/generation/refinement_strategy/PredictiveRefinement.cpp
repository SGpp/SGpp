// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "PredictiveRefinement.hpp"


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void PredictiveRefinement::collectRefinablePoints(
      GridStorage* storage, RefinementFunctor* functor,
      size_t refinements_num, size_t* max_indices,
      RefinementFunctor::value_type* max_values) {
    }

  } /* namespace base */
} /* namespace SGPP */
