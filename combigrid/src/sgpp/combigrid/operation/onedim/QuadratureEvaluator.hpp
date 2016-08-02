/*
 * QuadratureEvaluator.hpp
 *
 *  Created on: Feb 23, 2016
 *      Author: liedtkjn
 */

#ifndef QUADRATUREEVALUATOR_HPP_
#define QUADRATUREEVALUATOR_HPP_

#include "AbstractLinearEvaluator.hpp"
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/SingleFunction.hpp>
#include <vector>
#include <functional>

namespace SGPP {
namespace combigrid {

class QuadratureEvaluator: public AbstractLinearEvaluator<FloatScalarVector> {

	std::vector<SGPP::float_t> xValues;
	std::vector<FloatScalarVector> weights;
	SGPP::combigrid::SingleFunction weight_function;
	bool normalizeWeights;

public:
	QuadratureEvaluator(SGPP::combigrid::SingleFunction weight_function = SGPP::combigrid::SingleFunction(constantFunction<float_t>(static_cast<float_t>(1.0))), bool normalizeWeights = false);
	QuadratureEvaluator(QuadratureEvaluator const &other);
	virtual ~QuadratureEvaluator();

	virtual std::vector<FloatScalarVector> getBasisCoefficients() {
		return weights;
	}

	virtual void setGridPoints(std::vector<SGPP::float_t> const &newXValues);

	virtual bool needsOrderedPoints();
	virtual bool needsParameter();
	virtual void setParameter(FloatScalarVector const &param);

	virtual std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear();
};

#endif /* QUADRATUREEVALUATOR_HPP_ */

} /* namespace combigrid */
} /* namespace SGPP */
