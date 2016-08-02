/*
 * CombigridOperation.hpp
 *
 *  Created on: 04.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_COMBIGRIDOPERATION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_COMBIGRIDOPERATION_HPP_

#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include "onedim/AbstractLinearEvaluator.hpp"
#include <sgpp/combigrid/MultiFunction.hpp>

#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <cstddef>
#include <vector>
#include <memory>

namespace SGPP {
namespace combigrid {

class CombigridOperationImpl; // we use pimpl for not having to include all the template stuff in the header

class CombigridOperation {
	std::shared_ptr<CombigridOperationImpl> impl; // unique_ptr would be possible, but gives SWIG errors

public:
	CombigridOperation(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
			MultiFunction func);

	CombigridOperation(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
			std::shared_ptr<AbstractCombigridStorage> storage);

	// TODO: add extra functions, for example for configuring the storage

	SGPP::float_t evaluate(size_t q, base::DataVector const &param = base::DataVector(0));

	// TODO: add static constructor functions
	static std::shared_ptr<CombigridOperation> createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
			MultiFunction func);
	static std::shared_ptr<CombigridOperation> createExpLejaPolynomialInterpolation(size_t numDimensions,
				MultiFunction func);
	static std::shared_ptr<CombigridOperation> createExpUniformPolynomialInterpolation(size_t numDimensions,
					MultiFunction func);
	static std::shared_ptr<CombigridOperation> createLinearClenshawCurtisPolynomialInterpolation(size_t numDimensions,
				MultiFunction func);
	static std::shared_ptr<CombigridOperation> createLinearLejaPolynomialInterpolation(size_t numDimensions,
				MultiFunction func);
	static std::shared_ptr<CombigridOperation> createLinearUniformPolynomialInterpolation(size_t numDimensions,
						MultiFunction func);
	static std::shared_ptr<CombigridOperation> createLinearLejaQuadrature(size_t numDimensions,
						MultiFunction func, size_t growthFactor = 2);
		static std::shared_ptr<CombigridOperation> createExpUniformLinearInterpolation(size_t numDimensions,
						MultiFunction func);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_COMBIGRIDOPERATION_HPP_ */
