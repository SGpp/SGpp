/*
 * AbstractFullGridEvaluator.hpp
 *
 *  Created on: 11.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTFULLGRIDEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTFULLGRIDEVALUATOR_HPP_

#include "../../definitions.hpp"
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>

#include <memory>

namespace sgpp{
namespace combigrid {

template<typename V> class AbstractFullGridEvaluator {
public:
	virtual ~AbstractFullGridEvaluator(){}

	virtual V eval(MultiIndex const &level) = 0;

	virtual size_t maxNewPoints(MultiIndex const &level) = 0;

	virtual std::shared_ptr<AbstractCombigridStorage> getStorage() = 0;

	virtual std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level, ThreadPool::Task callback, std::mutex &lock) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTFULLGRIDEVALUATOR_HPP_ */
