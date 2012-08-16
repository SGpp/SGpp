/*
 * AdaptiveSerialCombiGrid.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: kowitz_local
 */

#include "combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"

std::vector<int> combigrid::AdaptiveSerialCombiGrid::addToCombiScheme(
		std::vector<int> level) {

	COMBIGRID_ERROR_TEST((int)level.size()==abs(combikernel_->getDim()),
			"new level vector has the wrong size!")
	CombigridLevelVector current(combischeme_->getLevels(),
			combischeme_->getCoef());
	int currentSize = combischeme_->getNrSapces();
//	current.printLevelVec();
//	std::cout<<"-----------"<<std::endl;
	current = current.getChanges(level);
//	current.printLevelVec();
	std::vector<int> changes = combischeme_->updateScheme(current.getLevelVec(),
			current.getCoef());
	std::vector<int> newCoef(changes.size());

	std::vector<std::vector<int> > newLevels = current.getLevelVec();
	for (unsigned int i = 0; i < changes.size(); ++i) {
		newCoef[i] = combischeme_->getCoef(changes[i]);
	}
	(*combikernel_).updateCombiScheme(/*combischeme_*/newCoef, newLevels,
			changes);

//	for (int i = 0; i < newLevels.size(); i++) {
//		for (int j = 0; j < newLevels[i].size(); j++)
//			std::cout << newLevels[i][j] << '\t';
//		std::cout << std::endl;
//	}
//	std::cout << "changes: ";
//	for (int i = 0; i < changes.size(); i++)
//		std::cout << changes[i] << '\t';
//	std::cout << std::endl;

	std::vector<int> newGridsIndices;
	for (unsigned int i = 0; i < changes.size(); i++) {
		if (changes[i] >= currentSize)
			newGridsIndices.push_back(changes[i]);
	}
	return newGridsIndices;

//	for ( int i = 0; i < combischeme_->getNrSapces(); i++) {
//		for ( int j = 0; j < combischeme_->getDim(); j++) {
//			std::cout << combischeme_->getLevel(i)[j] << "\t";
//		}
//		std::cout << std::endl;
//	}
//	for ( int i = 0; i < combischeme_->getNrSapces(); i++) {
//		for ( int j = 0; j < combischeme_->getDim(); j++) {
//			std::cout << combikernel_->getFullGridLevel(i)[j] << "\t";
//		}
//		std::cout << std::endl;
//	}

}
