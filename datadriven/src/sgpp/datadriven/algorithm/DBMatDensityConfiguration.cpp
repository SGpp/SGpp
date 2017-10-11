// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

namespace sgpp {
namespace datadriven {

DBMatDensityConfiguration::DBMatDensityConfiguration()
    : grid_type_(sgpp::base::GridType::Linear),
      grid_dim_(0),
      grid_level_(0),
      numRefinements_(0),
      ref_threshold_(0),
      ref_noPoints_(0),
      regularization_(RegularizationType::Identity),
      lambda_(0),
      decomp_type_(DBMatDecompostionType::Chol) {}

DBMatDensityConfiguration::DBMatDensityConfiguration(const sgpp::base::RegularGridConfiguration& gc,
                                                     const sgpp::base::AdpativityConfiguration& ac,
                                                     sgpp::datadriven::RegularizationType reg,
                                                     double lambda, DBMatDecompostionType dt)
    : grid_type_(gc.type_),
      grid_dim_(gc.dim_),
      grid_level_(gc.level_),
      numRefinements_(ac.numRefinements_),
      ref_threshold_(ac.threshold_),
      ref_noPoints_(ac.noPoints_),
      regularization_(reg),
      lambda_(lambda),
      decomp_type_(dt) {}

}  // namespace datadriven
}  // namespace sgpp
