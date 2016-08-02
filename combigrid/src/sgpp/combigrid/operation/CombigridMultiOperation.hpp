/*
 * CombigridMultiOperation.hpp
 *
 *  Created on: 06.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_COMBIGRIDMULTIOPERATION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_COMBIGRIDMULTIOPERATION_HPP_
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/algebraic/ArrayVector.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include "onedim/AbstractLinearEvaluator.hpp"
#include <sgpp/combigrid/MultiFunction.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <cstddef>
#include <vector>
#include <memory>

namespace sgpp{
namespace combigrid {

class CombigridMultiOperationImpl; // we use pimpl for not having to include all the template stuff in the header

class CombigridMultiOperation {
	std::shared_ptr<CombigridMultiOperationImpl> impl; // unique_ptr causes SWIG errors

public:
	CombigridMultiOperation(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
			MultiFunction func);

	CombigridMultiOperation(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
			std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
			std::shared_ptr<AbstractCombigridStorage> storage);

	// TODO: add extra functions, for example for configuring the storage

	base::DataVector evaluate(size_t q, std::vector<base::DataVector> const &params); // TODO: maybe change to DataMatrix?
	base::DataVector evaluateAdaptive(size_t maxNumPoints, std::vector<base::DataVector> const &params); // TODO: maybe change to DataMatrix?

	std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> getDifferences();

	// TODO: add static constructor functions
	static std::shared_ptr<CombigridMultiOperation> createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
			MultiFunction func);
	static std::shared_ptr<CombigridMultiOperation> createLinearClenshawCurtisPolynomialInterpolation(size_t numDimensions,
			MultiFunction func);
	static std::shared_ptr<CombigridMultiOperation> createExpLejaPolynomialInterpolation(size_t numDimensions,
					MultiFunction func);
	static std::shared_ptr<CombigridMultiOperation> createLinearLejaPolynomialInterpolation(size_t numDimensions,
				MultiFunction func, size_t growthFactor = 2);
	static std::shared_ptr<CombigridMultiOperation> createLinearLejaQuadrature(size_t numDimensions,
					MultiFunction func, size_t growthFactor = 2);
	static std::shared_ptr<CombigridMultiOperation> createExpUniformLinearInterpolation(size_t numDimensions,
					MultiFunction func);
	// static std::shared_ptr<CombigridOperation> createLejaLinearPolynomialInterpolation(size_t numDimensions); // TODO
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_COMBIGRIDMULTIOPERATION_HPP_ */
