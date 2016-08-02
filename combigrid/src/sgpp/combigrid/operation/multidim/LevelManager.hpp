/*
 * LevelManager.hpp
 *
 *  Created on: 22.07.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELMANAGER_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELMANAGER_HPP_

#include "../../definitions.hpp"
#include "../../common/BoundedSumMultiIndexIterator.hpp"
#include "../../storage/AbstractMultiStorage.hpp"
#include "../../storage/tree/TreeStorage.hpp"
#include "AdaptiveRefinementStrategy.hpp"
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include "LevelHelpers.hpp"
#include <memory>
#include <queue>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <mutex>

#include "CombigridEvaluator.hpp" // TODO: remove
#include <sgpp/combigrid/algebraic/ScalarVector.hpp> // TODO: remove

namespace sgpp{
namespace combigrid {

// TODO: put functions in cpp file

class LevelManager {
protected:
	// data structures for adaptive refinement
	MultiIndexQueue queue;
	std::shared_ptr<TreeStorage<std::shared_ptr<LevelInfo>>> levelData;
	size_t numDimensions;
	std::shared_ptr<CombigridEvaluator<FloatScalarVector>> combiEval;
	std::mutex managerMutex;

	/**
	 * By implementing this method in a derived class, the adaption can be customized.
	 */
	virtual double computePriority(MultiIndex const &level) = 0;

	/**
	 * Initializes the data structures for adaptive level generation
	 */
	virtual void initAdaption();

	virtual void tryAddSuccessors(MultiIndex const &level);

	virtual void tryAddLevel(MultiIndex const &level);

	virtual void addToQueue(MultiIndex const &level, std::shared_ptr<LevelInfo> levelInfo);

	virtual std::vector<MultiIndex> getPredecessors(MultiIndex const &level);

	virtual std::vector<MultiIndex> getSuccessors(MultiIndex const &level);

	virtual void beforeComputation(MultiIndex const &level);

	virtual void afterComputation(MultiIndex const &level);

	virtual void predecessorsCompleted(MultiIndex const &level);

	virtual void updatePriority(MultiIndex const &level, std::shared_ptr<LevelInfo> levelInfo);

	//---------- old methods ----------

	// If a level has been computed, we try to add its successors to the queue if all of their predecessors are computed.
	virtual void tryAddSuccessorsToQueue(MultiIndex const &index);

	// Adds a level to the queue if it is not in the CombigridEvaluator, but all of its predecessors are.
	virtual void tryAddIndexToQueue(MultiIndex const &nextIndex);

	/**
	 * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as in most papers).
	 * If you have a norm w with levels starting from 1, simply use q = w - dim().
	 */
	std::vector<MultiIndex> getRegularLevels(size_t q);

	std::vector<MultiIndex> getRegularLevelsByNumPoints(size_t maxNumPoints);

	void precomputeLevelsParallel(std::vector<MultiIndex> const &levels, size_t numThreads);

	void addLevels(std::vector<MultiIndex> const &levels);

public:
	// TODO constructor

	virtual ~LevelManager();

	/**
	 * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as in most papers).
	 * If you have a norm w with levels starting from 1, simply use q = w - dim().
	 */
	void addRegularLevels(size_t q);

	void addRegularLevelsByNumPoints(size_t maxNumPoints);

	/**
	 * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as in most papers).
	 * If you have a norm w with levels starting from 1, simply use q = w - dim().
	 */
	void addRegularLevelsParallel(size_t q, size_t numThreads);

	void addRegularLevelsByNumPointsParallel(size_t maxNumPoints, size_t numThreads);

	virtual size_t dim() const;

	std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() const;

	std::string getSerializedLevelStructure() const;

	void addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage);

	void addLevelsFromSerializedStructure(std::string serializedStructure);

	virtual void addLevelsAdaptive(size_t maxNumPoints);

	virtual void addLevelsAdaptiveParallel(size_t maxNumPoints, size_t numThreads);
};

}
/* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELMANAGER_HPP_ */
