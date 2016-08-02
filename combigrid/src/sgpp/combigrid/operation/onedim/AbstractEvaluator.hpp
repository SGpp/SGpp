/*
 * AbstractEvaluator.hpp
 *
 *  Created on: 11.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTEVALUATOR_HPP_

#include <sgpp/globaldef.hpp>
#include <vector>
#include <memory>

namespace SGPP {
namespace combigrid {

template<typename V> class AbstractEvaluator {
public:
	virtual ~AbstractEvaluator(){}

	virtual void setGridPoints(std::vector<SGPP::float_t> const &xValues) = 0;
	virtual std::shared_ptr<AbstractEvaluator<V>> clone() = 0;
	virtual bool needsOrderedPoints() = 0;
	virtual bool needsParameter() = 0;
	virtual void setParameter(V const &param) = 0;

	virtual V eval(std::vector<float_t> const &functionValues) = 0;

	virtual V eval(std::vector<V> const &functionValues) = 0;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTEVALUATOR_HPP_ */
