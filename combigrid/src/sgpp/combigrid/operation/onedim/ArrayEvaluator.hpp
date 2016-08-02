/*
 * ArrayEvaluator.hpp
 *
 *  Created on: 04.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_

#include "AbstractLinearEvaluator.hpp"
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/algebraic/ArrayVector.hpp>
#include <sgpp/combigrid/definitions.hpp>

namespace SGPP {
namespace combigrid {

template<typename ScalarEvaluator> class ArrayEvaluator: public AbstractLinearEvaluator<FloatArrayVector> {
	std::vector<ScalarEvaluator> evaluators;
	std::vector<FloatArrayVector> basisCoefficients;
	std::vector<SGPP::float_t> xValues;
	bool doesNeedParameter;

	void computeBasisCoefficients() {
		basisCoefficients = std::vector<FloatArrayVector>(xValues.size(), FloatArrayVector::zero());

		for(size_t i = 0; i < evaluators.size(); ++i) {
			auto coeff = evaluators[i].getBasisCoefficients();
			for(size_t j = 0; j < coeff.size(); ++j) {
				basisCoefficients[j].at(i) = coeff[j];
			}
		}
	}

public:
	ArrayEvaluator(bool doesNeedParameter, ScalarEvaluator evaluatorPrototype = ScalarEvaluator())
		: evaluators(1, evaluatorPrototype)
		, basisCoefficients()
		, xValues()
		, doesNeedParameter(doesNeedParameter) {
	}

	ArrayEvaluator(ArrayEvaluator<ScalarEvaluator> const &other)
		: evaluators(other.evaluators)
		, basisCoefficients(other.basisCoefficients)
		, xValues(other.xValues)
		, doesNeedParameter(other.doesNeedParameter) {

	}

	~ArrayEvaluator() {

	}

	virtual std::vector<FloatArrayVector> getBasisCoefficients() {
		return basisCoefficients;
	}

	virtual void setGridPoints(std::vector<SGPP::float_t> const &xValues) {
		this->xValues = xValues;
		for(auto &eval : evaluators) {
			eval.setGridPoints(xValues);
		}

		computeBasisCoefficients();
	}

	// don't change return type back to the equivalent std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> or swig will give compile errors...
	virtual std::shared_ptr<AbstractLinearEvaluator<ArrayVector<float_t, ScalarVector<float_t>>>> cloneLinear() {
		return std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>(new ArrayEvaluator<ScalarEvaluator>(*this));
	}

	virtual bool needsOrderedPoints() {
		return evaluators[0].needsOrderedPoints();
	}

	virtual bool needsParameter() {
		return doesNeedParameter;
	}

	virtual void setParameter(FloatArrayVector const &param) {
		if(evaluators.size() > param.size()) {
			evaluators.resize(param.size());
		} else {
			if(evaluators.size() == 0) {
				evaluators.emplace_back();
				evaluators.back().setGridPoints(xValues);
			}
			while(evaluators.size() < param.size()) {
				evaluators.push_back(evaluators.back());
			}
		}
		for(size_t i = 0; i < param.size(); ++i) {
			evaluators[i].setParameter(param[i]);
		}

		computeBasisCoefficients();
	}
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_ */
