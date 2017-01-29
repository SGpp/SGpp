// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

namespace sgpp {
namespace datadriven {

DBMatDensityConfiguration::DBMatDensityConfiguration(
    sgpp::base::RegularGridConfiguration* gc,
    sgpp::base::AdpativityConfiguration* ac,
    sgpp::datadriven::RegularizationType reg, double lambda,
    DBMatDecompostionType dt)
    : grid_type_(gc->type_),
      grid_dim_(gc->dim_),
      grid_level_(gc->level_),
      numRefinements_(ac->numRefinements_),
      ref_threshold_(ac->threshold_),
      ref_noPoints_(ac->noPoints_),
      regularization_(reg),
      lambda_(lambda),
      decomp_type_(dt) {}

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
