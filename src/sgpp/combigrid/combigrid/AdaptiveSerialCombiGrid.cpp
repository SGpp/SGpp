/*
 * AdaptiveSerialCombiGrid.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: kowitz_local
 */

#include "combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"


void combigrid::AdaptiveSerialCombiGrid::addToCombiScheme(
		std::vector<int> level) {
	CombigridLevelVector current(combischeme_->getLevels(),
			combischeme_->getCoef());
	current = current.getChanges(level);
	std::vector<int> changes = combischeme_->updateScheme(
			current.getLevelVec(), current.getCoef());
	std::vector<int> newCoef(changes.size());
	for (unsigned int i = 0; i < changes.size(); ++i) {
		newCoef[i] = combischeme_->getCoef(changes[i]);
	}
	(*combikernel_).updateCombiScheme(/*combischeme_*/newCoef, level, changes);
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
