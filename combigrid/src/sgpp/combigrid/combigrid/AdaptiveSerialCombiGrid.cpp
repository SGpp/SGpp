// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/combigrid/combigrid/AdaptiveSerialCombiGrid.hpp>

std::vector<int> combigrid::AdaptiveSerialCombiGrid::addToCombiScheme(
  std::vector<int> level) {

  COMBIGRID_ERROR_TEST((int)level.size() == abs(combikernel_->getDim()),
                       "new level vector has the wrong size!")
  CombigridLevelVector current(combischeme_->getLevels(),
                               combischeme_->getCoef());
  int currentSize = combischeme_->getNrSapces();
  current = current.getChanges(level);
  std::vector<int> changes = combischeme_->updateScheme(current.getLevelVec(),
                             current.getCoef());
  std::vector<double> newCoef(changes.size());
  std::vector<std::vector<int> > newLevels = current.getLevelVec();

  for (unsigned int i = 0; i < changes.size(); ++i) {
    newCoef[i] = combischeme_->getCoef(changes[i]);
  }

  (*combikernel_).updateCombiScheme(/*combischeme_*/newCoef, newLevels,
      changes);


  std::vector<int> newGridsIndices;

  for (unsigned int i = 0; i < changes.size(); i++) {
    if (changes[i] >= currentSize)
      newGridsIndices.push_back(changes[i]);
  }

  return newGridsIndices;


}