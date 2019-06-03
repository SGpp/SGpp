// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "LevelManager.hpp"

#include <sgpp/combigrid/threading/PtrGuard.hpp>

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

LevelManager::LevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator,
                           bool collectStats)
    : queue(),
      numDimensions(levelEvaluator->numDims()),
      combiEval(levelEvaluator),
      collectStats(collectStats) {
  managerMutex = std::make_shared<std::recursive_mutex>();
  infoOnAddedLevels = std::make_shared<LevelInfos>();
  levelData = std::make_shared<TreeStorage<std::shared_ptr<LevelInfo>>>(numDimensions);
}

LevelManager::LevelManager() : queue(), numDimensions(0), combiEval(nullptr), collectStats(false) {
  managerMutex = std::make_shared<std::recursive_mutex>();
  infoOnAddedLevels = std::make_shared<LevelInfos>();
  levelData = std::make_shared<TreeStorage<std::shared_ptr<LevelInfo>>>(numDimensions);
}

LevelManager::~LevelManager() {}

void LevelManager::setLevelEvaluator(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator) {
  combiEval = levelEvaluator;
  numDimensions = levelEvaluator->numDims();
}

void LevelManager::initAdaption() {
  queue.clear();
  levelData.reset(new TreeStorage<std::shared_ptr<LevelInfo>>(numDimensions));

  auto levelStructure = getLevelStructure();
  auto it = levelStructure->getStoredDataIterator();

  while (it->isValid()) {
    auto level = it->getMultiIndex();
    levelData->set(level, std::make_shared<LevelInfo>(combiEval->getDifferenceNorm(level),
                                                      combiEval->maxNewPoints(level),
                                                      combiEval->numPoints(level)));
    it->moveToNext();
  }

  auto dataIt = levelData->getStoredDataIterator();
  while (dataIt->isValid()) {
    tryAddSuccessors(dataIt->getMultiIndex());
    dataIt->moveToNext();
  }

  if (queue.empty()) {
    // no difference has yet been computed, so add start value
    MultiIndex startLevel(numDimensions, 0);
    tryAddLevel(startLevel);
  }
}

void LevelManager::tryAddSuccessors(const MultiIndex &level) {
  for (auto succLevel : getSuccessors(level)) {
    tryAddLevel(succLevel);
  }
}

void LevelManager::tryAddLevel(const MultiIndex &level) {
  if (levelData->containsIndex(level)) {
    return;
  }

  auto predecessors = getPredecessors(level);
  size_t numPredecessors = predecessors.size();
  size_t numMissingPredecessors = numPredecessors;

  for (auto &predLevel : predecessors) {
    if (levelData->containsIndex(predLevel)) {
      auto levelInfo = levelData->get(predLevel);
      if (levelInfo->computationStage >= ComputationStage::STARTED) {
        --numMissingPredecessors;
      }
    }
  }

  if (numMissingPredecessors > 0) {
    return;
  }

  auto levelInfo = std::make_shared<LevelInfo>(numPredecessors, combiEval->maxNewPoints(level),
                                               combiEval->numPoints(level));

  for (auto &predLevel : predecessors) {
    auto predInfo = levelData->get(predLevel);
    if (predInfo->computationStage == ComputationStage::COMPLETED) {
      --levelInfo->numNotStartedPredecessors;
      --levelInfo->numNotCompletedPredecessors;
    } else if (predInfo->computationStage == ComputationStage::STARTED ||
               predInfo->computationStage == ComputationStage::TERMINATED) {
      --levelInfo->numNotStartedPredecessors;
    }
  }

  levelData->set(level, levelInfo);

  if (levelInfo->numNotStartedPredecessors == 0) {
    // the new level may be started, so we can add it to the priority queue
    addToQueue(level, levelInfo);
  }
}

void LevelManager::addToQueue(const MultiIndex &level, std::shared_ptr<LevelInfo> levelInfo) {
  double priority = computePriority(level);
  size_t maxNewPoints = combiEval->maxNewPoints(level);
  QueueEntry entry(level, priority, maxNewPoints);
  auto handle = queue.push(entry);
  levelInfo->handle = std::make_shared<MultiIndexQueue::handle_type>(handle);
}

void LevelManager::beforeComputation(const MultiIndex &level) {
  auto levelInfo = levelData->get(level);
  levelInfo->computationStage = ComputationStage::STARTED;
  /*
   * Invalidate the handle since the element has been popped from the queue.
   */
  levelInfo->handle = nullptr;

  auto successors = getSuccessors(level);
  for (auto &succLevel : successors) {
    if (levelData->containsIndex(succLevel)) {
      auto succInfo = levelData->get(succLevel);
      --succInfo->numNotStartedPredecessors;

      if (succInfo->numNotStartedPredecessors == 0) {
        // the new level may be started, so we can add it to the priority queue
        addToQueue(succLevel, succInfo);
      }
    } else {
      tryAddLevel(succLevel);
    }
  }
}

void LevelManager::afterComputation(const MultiIndex &level) {
  auto levelInfo = levelData->get(level);
  levelInfo->computationStage = ComputationStage::TERMINATED;

  if (levelInfo->numNotCompletedPredecessors == 0) {
    predecessorsCompleted(level);
  }
}

void LevelManager::predecessorsCompleted(const MultiIndex &level) {
  auto levelInfo = levelData->get(level);
  levelInfo->computationStage = ComputationStage::COMPLETED;

  bool validResult = addLevelToCombiEval(level);
  levelInfo->norm =
      validResult ? combiEval->getDifferenceNorm(level) : 0.0;       // TODO(holzmudd): improve?
  levelInfo->priority = validResult ? computePriority(level) : 0.0;  // TODO(franzefn): improve?
  auto successors = getSuccessors(level);
  for (auto &succLevel : successors) {
    if (!levelData->containsIndex(succLevel)) {
      continue;
    }
    auto succInfo = levelData->get(succLevel);
    --succInfo->numNotCompletedPredecessors;
    if (succInfo->computationStage == ComputationStage::TERMINATED &&
        succInfo->numNotCompletedPredecessors == 0) {
      predecessorsCompleted(succLevel);
    } else if (succInfo->handle != nullptr) {
      // If handle is nullptr, this means that the element has already been started and thus removed
      // from the queue.
      updatePriority(succLevel, succInfo);
    }
  }
}

void LevelManager::updatePriority(const MultiIndex &level, std::shared_ptr<LevelInfo> levelInfo) {
  double priority = computePriority(level);
  levelInfo->setPriority(queue, priority);
}

std::vector<MultiIndex> LevelManager::getPredecessors(MultiIndex const &level) {
  std::vector<MultiIndex> result;

  for (size_t prevDim = 0; prevDim < numDimensions; ++prevDim) {
    if (level[prevDim] > 0) {
      MultiIndex prevIndex = level;
      --prevIndex[prevDim];
      result.push_back(prevIndex);
    }
  }

  return result;
}

std::vector<MultiIndex> LevelManager::getSuccessors(MultiIndex const &level) {
  std::vector<MultiIndex> result;

  for (size_t dim = 0; dim < numDimensions; ++dim) {
    MultiIndex nextIndex = level;
    ++nextIndex[dim];
    result.push_back(nextIndex);
  }

  return result;
}

std::vector<MultiIndex> LevelManager::getRegularLevels(size_t q) {
  std::vector<MultiIndex> result;

  BoundedSumMultiIndexIterator iterator(numDimensions, q);

  while (iterator.isValid()) {
    result.push_back(iterator.value());
    iterator.moveToNext();
  }

  return result;
}

std::vector<MultiIndex> LevelManager::getRegularLevelsByNumPoints(size_t maxNumPoints) {
  std::vector<MultiIndex> result;

  size_t numPoints = 0;

  size_t q = 0;

  while (true) {
    BoundedSumMultiIndexIterator iterator(numDimensions, q);

    while (iterator.isValid()) {
      MultiIndex nextLevel = iterator.value();

      if (!combiEval->containsLevel(nextLevel)) {
        size_t maxNewPoints = combiEval->maxNewPoints(nextLevel);
        numPoints += maxNewPoints;

        if (numPoints > maxNumPoints) {
          return result;
        }

        result.push_back(nextLevel);
      }

      iterator.moveToNext();
    }

    ++q;
  }
}

void LevelManager::precomputeLevelsParallel(const std::vector<MultiIndex> &levels,
                                            size_t numThreads) {
  auto threadPool = std::make_shared<ThreadPool>(numThreads, ThreadPool::terminateWhenIdle);
  combiEval->setMutex(managerMutex);
  for (auto &level : levels) {
    threadPool->addTasks(combiEval->getLevelTasks(level, ThreadPool::Task([]() {})));
  }
  threadPool->start();
  threadPool->join();
  combiEval->setMutex(nullptr);
}

void LevelManager::addStats(const MultiIndex &level) {
  // load level info. If not existing, load it
  LevelInfo levelInfo(combiEval->getDifferenceNorm(level), combiEval->maxNewPoints(level),
                      combiEval->numPoints(level));
  infoOnAddedLevels->insert(level, levelInfo);
}

bool LevelManager::addLevelToCombiEval(const MultiIndex &level) {
  bool isOldSubspace = combiEval->containsLevel(level);
  bool isValid = combiEval->addLevel(level);
  if (collectStats && !isOldSubspace) {
    addStats(level);
  }
  return isValid;
}

void LevelManager::addLevels(const std::vector<MultiIndex> &levels) {
  for (MultiIndex level : levels) {
    addLevelToCombiEval(level);
  }
}

void LevelManager::addRegularLevelsParallel(size_t q, size_t numThreads) {
  auto levels = getRegularLevels(q);
  precomputeLevelsParallel(levels, numThreads);
  // update stats vector
  infoOnAddedLevels->incrementCounter();
  addLevels(levels);
}

void LevelManager::addRegularLevelsByNumPointsParallel(size_t maxNumPoints, size_t numThreads) {
  auto levels = getRegularLevelsByNumPoints(maxNumPoints);
  precomputeLevelsParallel(levels, numThreads);
  // update stats
  infoOnAddedLevels->incrementCounter();
  addLevels(levels);
}

void LevelManager::addRegularLevels(size_t q) {
  auto levels = getRegularLevels(q);
  // update stats vector
  infoOnAddedLevels->incrementCounter();
  addLevels(levels);
}

void LevelManager::addRegularLevelsByNumPoints(size_t maxNumPoints) {
  auto levels = getRegularLevelsByNumPoints(maxNumPoints);
  // update stats vector
  infoOnAddedLevels->incrementCounter();
  addLevels(levels);
}

size_t LevelManager::numDims() const { return numDimensions; }

std::shared_ptr<TreeStorage<uint8_t>> LevelManager::getLevelStructure() const {
  return combiEval->getLevelStructure();
}

std::string LevelManager::getSerializedLevelStructure() const {
  return TreeStorageSerializationStrategy<uint8_t>(numDimensions).serialize(getLevelStructure());
}

void LevelManager::addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage) {
  if (storage != nullptr) {
    auto it = storage->getStoredDataIterator();
    // update stats vector
    infoOnAddedLevels->incrementCounter();
    while (it->isValid()) {
      addLevelToCombiEval(it->getMultiIndex());
      it->moveToNext();
    }
  }
}

void LevelManager::addLevelsFromStructureParallel(std::shared_ptr<TreeStorage<uint8_t>> storage,
                                                  size_t numThreads) {
  if (storage != nullptr) {
    auto it = storage->getStoredDataIterator();

    std::vector<MultiIndex> levels;
    while (it->isValid()) {
      levels.push_back(it->getMultiIndex());
      it->moveToNext();
    }
    precomputeLevelsParallel(levels, numThreads);
    infoOnAddedLevels->incrementCounter();
    addLevels(levels);
  }
}

void LevelManager::addLevelsFromSerializedStructure(std::string serializedStructure) {
  addLevelsFromStructure(
      TreeStorageSerializationStrategy<uint8_t>(numDimensions).deserialize(serializedStructure));
}

void LevelManager::addLevelsFromSerializedStructureParallel(std::string serializedStructure,
                                                            size_t numThreads) {
  addLevelsFromStructureParallel(
      TreeStorageSerializationStrategy<uint8_t>(numDimensions).deserialize(serializedStructure),
      numThreads);
}

void LevelManager::addLevelsAdaptive(size_t maxNumPoints, bool verbose) {
  initAdaption();

  size_t currentPointBound = 0;

  infoOnAddedLevels->incrementCounter();
  while (!queue.empty()) {
    // print current queue
    //    queue.print();

    QueueEntry entry = queue.top();
    queue.pop();

    currentPointBound += entry.maxNewPoints;

    if (currentPointBound > maxNumPoints) {
      break;
    }

    beforeComputation(entry.level);  // successors may be added here
    afterComputation(entry.level);   // the actual computation is done here, in addLevel() and the
                                     // successors are updated here
    if (verbose) {
      std::cout << "added level ";
      for (auto &index : entry.level) {
        std::cout << index << " ";
      }
      std::cout << "\n";
    }
  }
}

void LevelManager::addLevelsAdaptiveParallel(size_t maxNumPoints, size_t numThreads) {
  initAdaption();

  size_t currentPointBound = 0;
  infoOnAddedLevels->incrementCounter();

  combiEval->setMutex(managerMutex);

  auto threadPool = std::make_shared<ThreadPool>(
      numThreads,
      ThreadPool::IdleCallback([&currentPointBound, maxNumPoints, this](ThreadPool &tp) {
        CGLOG_SURROUND(PtrGuard guard(managerMutex));
        if (queue.empty()) {
          std::cout << "Error: queue is empty\n";
          CGLOG("leave guard(*managerMutex)");
          return;
        }

        // print current queue
        //        queue.print();

        QueueEntry entry = queue.top();

        currentPointBound += entry.maxNewPoints;

        if (currentPointBound > maxNumPoints) {
          tp.triggerTermination();
          CGLOG("leave guard(*managerMutex)");
          return;
        }

        CGLOG("before beforeComputation()");
        /*
         * Must be placed after the triggerTermination() in the if clause since after the pop()
         * operation, the corresponding handle in the LevelInfo must be set to nullptr in order to
         * avoid accessing a handle to a popped element. Invalidating the handle is done by
         * beforeComputation().
         */
        queue.pop();
        beforeComputation(entry.level);
        CGLOG("before getLevelTasks()");
        auto tasks = combiEval->getLevelTasks(entry.level, ThreadPool::Task([this, entry]() {
                                                // the mutex will be locked when this callback is
                                                // called
                                                // PtrGuard guard(this->managerMutex);
                                                afterComputation(entry.level);
                                              }));
        CGLOG("before addTasks()");
        tp.addTasks(tasks);
        CGLOG("leave guard(*managerMutex)");
      }));

  threadPool->start();
  threadPool->join();

  combiEval->setMutex(nullptr);
}

void LevelManager::addLevelsAdaptiveByNumLevels(size_t numLevels) {
  initAdaption();

  infoOnAddedLevels->incrementCounter();
  for (size_t i = 0; i < numLevels; ++i) {
    QueueEntry entry = queue.top();
    queue.pop();

    beforeComputation(entry.level);  // successors may be added here
    afterComputation(entry.level);   // the actual computation is done here
  }
}

sgpp::base::DataMatrix LevelManager::convertLevelStructureToMatrix(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const &levelstructure, size_t numDims) {
  sgpp::base::DataMatrix levelstructureMatrix(0, numDims);
  auto it = levelstructure->getStoredDataIterator();
  while (it->isValid()) {
    sgpp::combigrid::MultiIndex index = it->getMultiIndex();
    sgpp::base::DataVector row;
    for (auto &i : index) {
      row.push_back(static_cast<double>(i));
    }
    levelstructureMatrix.appendRow(row);
    it->moveToNext();
  }
  return levelstructureMatrix;
}

void LevelManager::printLevelStructure(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const &levelstructure) {
  auto it = levelstructure->getStoredDataIterator();
  while (it->isValid()) {
    sgpp::combigrid::MultiIndex index = it->getMultiIndex();
    for (auto &i : index) {
      std::cout << i << " ";
    }
    std::cout << "\n";
    it->moveToNext();
  }
}

void LevelManager::enableStatsCollection() { collectStats = true; }

void LevelManager::disableStatsCollection() { collectStats = false; }

std::shared_ptr<LevelInfos> LevelManager::getInfoOnAddedLevels() { return infoOnAddedLevels; }

} /* namespace combigrid */
} /* namespace sgpp*/
