/*
 * LinearInterpolationEvaluator.hpp
 *
 *  Created on: 16.02.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_LINEARINTERPOLATIONEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_LINEARINTERPOLATIONEVALUATOR_HPP_

#include "AbstractLinearEvaluator.hpp"
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>

namespace SGPP {
namespace combigrid {

class LinearInterpolationEvaluator: public AbstractLinearEvaluator<FloatScalarVector> {
	SGPP::float_t evaluationPoint;
	std::vector<FloatScalarVector> basisCoefficients;
	std::vector<SGPP::float_t> xValues;

	void computeBasisCoefficients();

public:
	LinearInterpolationEvaluator();
	virtual ~LinearInterpolationEvaluator();
	LinearInterpolationEvaluator(LinearInterpolationEvaluator const &other);

	virtual std::vector<FloatScalarVector> getBasisCoefficients() {
		return basisCoefficients;
	}

	virtual void setGridPoints(std::vector<SGPP::float_t> const &newXValues);
	virtual std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear();
	virtual bool needsOrderedPoints();
	virtual bool needsParameter();
	virtual void setParameter(FloatScalarVector const &param);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_LINEARINTERPOLATIONEVALUATOR_HPP_ */
