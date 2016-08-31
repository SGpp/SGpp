// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DBMATDENSITYCONFIGURATION_H_
#define DBMATDENSITYCONFIGURATION_H_

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

enum DBMatDecompostionType
{
  DBMatDecompLU, DBMatDecompEigen, DBMatDecompChol
};

/**
 * Class that stores all the configuration information 
 * for an offline object for classification with the 
 * density based approach
 */

namespace sgpp {
namespace datadriven {

class DBMatDensityConfiguration
{
public:
        
        /**
	 * Constructor for hierarchical basis grids
	 *
	 * @param rg the configuration for the sparse grid
	 * @param ac the adaptivity configuration 
	 * @param reg the type of regularization
	 * @param lambda the weighting factor for the regularization term (can be changed later for certain decomposition types e.g. Eigendecompostion)
	 * @param decomp the kind of decomposition that should be used
	 */
	DBMatDensityConfiguration(sgpp::base::RegularGridConfiguration* gc,
			          sgpp::base::AdpativityConfiguration* ac,
			          sgpp::datadriven::RegularizationType reg, 
                                  double lambda,
			          DBMatDecompostionType decomp);

	sgpp::base::GridType grid_type_; // grid type
	size_t grid_dim_; // number of dimensions
	size_t grid_level_; // grid_level (only for hierarchical basis grids)

        //REFINEMENT
	size_t numRefinements_; // number of refinements
	double ref_threshold_; // refinement threshold for surpluses
	//bool ref_maxLevelType_; // refinement type: false: classic, true: maxLevel
	size_t ref_noPoints_; // max. number of points to be refined
	//double ref_percent_; // max. percent of points to be refined

	//REGULARIZATION:
	sgpp::datadriven::RegularizationType regularization_; //regularization operator
	double lambda_; //regularization parameter lambda

	//DECOMPOSITION:
	DBMatDecompostionType decomp_type_; //Type of matrix decomposition

};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DBMATDENSITYCONFIGURATION_H_ */

