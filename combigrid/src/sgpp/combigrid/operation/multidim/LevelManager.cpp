// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "LevelManager.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

LevelManager::LevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator)
    : queue(),
      levelData(),
      numDimensions(levelEvaluator->dim()),
      combiEval(levelEvaluator),
      managerMutex(new std::mutex()) {}

LevelManager::~LevelManager() {}

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

  bool validResult = combiEval->addLevel(level);
  levelInfo->norm =
      validResult ? combiEval->getDifferenceNorm(level) : 0.0;  // TODO(holzmudd): improve?
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

      size_t maxNewPoints = combiEval->maxNewPoints(nextLevel);
      numPoints += maxNewPoints;

      if (numPoints > maxNumPoints) {
        return result;
      }

      result.push_back(nextLevel);
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
    threadPool->addTasks(combiEval->getLevelTasks(level, []() {}));
  }
  threadPool->start();
  threadPool->join();
  combiEval->setMutex(nullptr);
}

void LevelManager::addLevels(const std::vector<MultiIndex> &levels) {
  for (auto &level : levels) {
    combiEval->addLevel(level);
  }
}

void LevelManager::addRegularLevelsParallel(size_t q, size_t numThreads) {
  auto levels = getRegularLevels(q);
  precomputeLevelsParallel(levels, numThreads);
  addLevels(levels);
}

void LevelManager::addRegularLevelsByNumPointsParallel(size_t maxNumPoints, size_t numThreads) {
  auto levels = getRegularLevelsByNumPoints(maxNumPoints);
  precomputeLevelsParallel(levels, numThreads);
  addLevels(levels);
}

void LevelManager::addRegularLevels(size_t q) {
  auto levels = getRegularLevels(q);
  addLevels(levels);
}

void LevelManager::addRegularLevelsByNumPoints(size_t maxNumPoints) {
  auto levels = getRegularLevelsByNumPoints(maxNumPoints);
  addLevels(levels);
}

size_t LevelManager::dim() const { return numDimensions; }

std::shared_ptr<TreeStorage<uint8_t>> LevelManager::getLevelStructure() const {
  return combiEval->getLevelStructure();
}

std::string LevelManager::getSerializedLevelStructure() const {
  return TreeStorageSerializationStrategy<uint8_t>(numDimensions).serialize(getLevelStructure());
}

void LevelManager::addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage) {
  auto it = storage->getStoredDataIterator();

  while (it->isValid()) {
    combiEval->addLevel(it->getMultiIndex());
    it->moveToNext();
  }
}

void LevelManager::addLevelsFromSerializedStructure(std::string serializedStructure) {
  addLevelsFromStructure(
      TreeStorageSerializationStrategy<uint8_t>(numDimensions).deserialize(serializedStructure));
}

void LevelManager::addLevelsAdaptive(size_t maxNumPoints) {
  initAdaption();

  size_t currentPointBound = 0;

  while (!queue.empty()) {
    QueueEntry entry = queue.top();
    queue.pop();

    currentPointBound += entry.maxNewPoints;

    if (currentPointBound > maxNumPoints) {
      break;
    }

    beforeComputation(entry.level);  // successors may be added here
    afterComputation(entry.level);   // the actual computation is done here, in addLevel() and the
                                     // successors are updated here
  }
}

void LevelManager::addLevelsAdaptiveParallel(size_t maxNumPoints, size_t numThreads) {
  initAdaption();

  size_t currentPointBound = 0;

  combiEval->setMutex(managerMutex);

  auto threadPool = std::make_shared<ThreadPool>(numThreads, [&](ThreadPool &tp) {
    CGLOG_SURROUND(std::lock_guard<std::mutex> guard(*managerMutex));
    if (queue.empty()) {
      std::cout << "Error: queue is empty\n";
      CGLOG("leave guard(*managerMutex)");
      return;
    }

    auto entry = queue.top();
    queue.pop();

    currentPointBound += entry.maxNewPoints;

    if (currentPointBound > maxNumPoints) {
      tp.triggerTermination();
      CGLOG("leave guard(*managerMutex)");
      return;
    }

    CGLOG("before beforeComputation()");

    beforeComputation(entry.level);
    CGLOG("before getLevelTasks()");
    auto tasks = combiEval->getLevelTasks(entry.level, [=]() {
      // the mutex will be locked when this callback is called
      afterComputation(entry.level);
    });
    CGLOG("before addTasks()");
    tp.addTasks(tasks);
    CGLOG("leave guard(*managerMutex)");
  });

  threadPool->start();
  threadPool->join();

  combiEval->setMutex(nullptr);
}

} /* namespace combigrid */
} /* namespace sgpp*/
