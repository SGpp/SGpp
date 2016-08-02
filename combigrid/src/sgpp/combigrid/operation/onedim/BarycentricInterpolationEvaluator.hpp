/*
 * BarycentricInterpolationEvaluator.hpp
 *
 *  Created on: 04.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BARYCENTRICINTERPOLATIONEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BARYCENTRICINTERPOLATIONEVALUATOR_HPP_

#include "AbstractLinearEvaluator.hpp"
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>

namespace SGPP {
namespace combigrid {

class BarycentricInterpolationEvaluator: public AbstractLinearEvaluator<FloatScalarVector> {
	SGPP::float_t evaluationPoint;
	std::vector<FloatScalarVector> basisCoefficients;
	std::vector<SGPP::float_t> wValues;
	std::vector<SGPP::float_t> xValues;

	void computeBasisCoefficients();

public:
	BarycentricInterpolationEvaluator();
	virtual ~BarycentricInterpolationEvaluator();
	BarycentricInterpolationEvaluator(BarycentricInterpolationEvaluator const &other);

	virtual std::vector<FloatScalarVector> getBasisCoefficients() {
		return basisCoefficients;
	}

	virtual void setGridPoints(std::vector<SGPP::float_t> const &newXValues);
	virtual std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear();
	virtual bool needsOrderedPoints();
	virtual bool needsParameter();
	virtual void setParameter(FloatScalarVector const &param);

	// TODO: eval could be optimized...
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BARYCENTRICINTERPOLATIONEVALUATOR_HPP_ */
