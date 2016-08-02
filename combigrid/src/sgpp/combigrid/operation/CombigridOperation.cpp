/*
 * CombigridOperation.cpp
 *
 *  Created on: 04.01.2016
 *      Author: david
 */

#include "CombigridOperation.hpp"
#include "multidim/CombigridEvaluator.hpp"
#include "multidim/FullGridTensorEvaluator.hpp"
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>

namespace SGPP {
namespace combigrid {

class CombigridOperationImpl {
public:
	CombigridOperationImpl(std::vector<std::shared_ptr<AbstractPointHierarchy> > pointHierarchies,
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > > evaluatorPrototypes,
			std::shared_ptr<AbstractCombigridStorage> storage) :
			storage(storage), fullGridEval(new FullGridTensorEvaluator<FloatScalarVector>(storage, evaluatorPrototypes, pointHierarchies)), combiEval(
					new CombigridEvaluator<FloatScalarVector>(pointHierarchies.size(), fullGridEval)) {

	}

	std::shared_ptr<AbstractCombigridStorage> storage;
	std::shared_ptr<FullGridTensorEvaluator<FloatScalarVector>> fullGridEval;
	std::shared_ptr<CombigridEvaluator<FloatScalarVector>> combiEval;
};

CombigridOperation::CombigridOperation(std::vector<std::shared_ptr<AbstractPointHierarchy> > pointHierarchies,
		std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > > evaluatorPrototypes,
		MultiFunction func) :
		impl(new CombigridOperationImpl(pointHierarchies, evaluatorPrototypes, std::shared_ptr<AbstractCombigridStorage>(new CombigridTreeStorage(pointHierarchies, func)))) {
}

CombigridOperation::CombigridOperation(std::vector<std::shared_ptr<AbstractPointHierarchy> > pointHierarchies,
		std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > > evaluatorPrototypes,
		std::shared_ptr<AbstractCombigridStorage> storage) :
		impl(new CombigridOperationImpl(pointHierarchies, evaluatorPrototypes, storage)) {
}

SGPP::float_t CombigridOperation::evaluate(size_t q, base::DataVector const &param) {
	std::vector<FloatScalarVector> scalars(param.getSize());
	for(size_t i = 0; i < param.getSize(); ++i) {
		scalars[i].value() = param[i];
	}

	impl->fullGridEval->setParameters(scalars);
	impl->combiEval->clear();
	impl->combiEval->addRegularLevels(q);

	return impl->combiEval->getValue().value();
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NestedPointHierarchy>(std::make_shared<ClenshawCurtisDistribution>(), std::make_shared<ExponentialLevelorderPointOrdering>())),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<BarycentricInterpolationEvaluator>()),
			func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpLejaPolynomialInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NestedPointHierarchy>(std::make_shared<LejaPointDistribution>(), std::make_shared<IdentityPointOrdering>(std::make_shared<ExponentialGrowthStrategy>(), false))),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<BarycentricInterpolationEvaluator>()),
			func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformPolynomialInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NestedPointHierarchy>(std::make_shared<UniformPointDistribution>(), std::make_shared<ExponentialLevelorderPointOrdering>())),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<BarycentricInterpolationEvaluator>()),
			func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearClenshawCurtisPolynomialInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NonNestedPointHierarchy>(std::make_shared<ClenshawCurtisDistribution>(), std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(2), true))),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<BarycentricInterpolationEvaluator>()),
			func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaPolynomialInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NestedPointHierarchy>(std::make_shared<LejaPointDistribution>(), std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(2), true))),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<BarycentricInterpolationEvaluator>()),
			func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearUniformPolynomialInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NonNestedPointHierarchy>(std::make_shared<UniformPointDistribution>(), std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(2), true))),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<BarycentricInterpolationEvaluator>()),
			func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaQuadrature(size_t numDimensions, MultiFunction func,
		size_t growthFactor) {
	return std::make_shared<CombigridOperation>(
					std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
						std::make_shared<NestedPointHierarchy>(std::make_shared<LejaPointDistribution>(), std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor), false))),
					std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
						std::make_shared<QuadratureEvaluator>()),
					func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformLinearInterpolation(size_t numDimensions,
		MultiFunction func) {
	return std::make_shared<CombigridOperation>(
			std::vector<std::shared_ptr<AbstractPointHierarchy> >(numDimensions,
				std::make_shared<NestedPointHierarchy>(std::make_shared<UniformPointDistribution>(), std::make_shared<ExponentialLevelorderPointOrdering>())),
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > >(numDimensions,
				std::make_shared<LinearInterpolationEvaluator>()),
			func);
}

} /* namespace combigrid */
} /* namespace SGPP */
