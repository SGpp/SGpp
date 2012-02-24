/*
 * CombiCombiSchemeBasis.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#include "combigrid/combischeme/CombiSchemeBasis.hpp"

using namespace std;

void combigrid::CombiSchemeBasis::removeDuplicates() {
	// this is Alize's code 1 to 1
	size_t i = 0;
	size_t j;
	std::vector<std::vector<int> >::iterator gt = levels_vector_.begin();
	std::vector<double>::iterator ct = cofficients_.begin();
	int size_fg = levels_vector_.size();
	while ((int) i < size_fg - 1) {
		j = i + 1;
		while ((int) j < size_fg) {
			if (levels_vector_[i] == levels_vector_[j]) {
				levels_vector_.erase(gt + j);
				cofficients_[i] += cofficients_[j];
				cofficients_.erase(ct + j);
				size_fg--;
			} else
				j++;
		}
		if (cofficients_[i] == 0.0) {
			levels_vector_.erase(gt + i);
			cofficients_.erase(ct + i);
			size_fg--;
			if (i == 0) {
				gt = levels_vector_.begin();
				ct = cofficients_.begin();
			}
		} else
			i++;
	}
}
std::vector<int> combigrid::CombiSchemeBasis::updateScheme(
		std::vector<std::vector<int> > levelsNew, std::vector<double> coef) {
	removeDuplicates();
	std::vector<int> result(0);
	for (unsigned int i = 0; i < levelsNew.size(); i++) {
		bool found = false;
		int found_ind = -1;
		for (unsigned int j = 0; j < levels_vector_.size(); j++) {
			bool equal = true;
			for (unsigned int k = 0; k < levels_vector_[j].size(); ++k) {
				if (levelsNew[i][k] != levels_vector_[j][k]) {
					equal = false;
					break;
				}

			}
			if (equal == true) {
				found = true;
				found_ind = j;
				result.push_back(j);
				break;
			}
		}
		if (found) {
			cofficients_[found_ind] += coef[i];
		} else {
			levels_vector_.push_back(levelsNew[i]);
			cofficients_.push_back(coef[i]);
			result.push_back(cofficients_.size() - 1);
		}
	}
	return result;
}

void combigrid::CombiSchemeBasis::setCoef(std::vector<double> newCoef) {
	COMBIGRID_ERROR_TEST((newCoef.size()==cofficients_.size()),
			"New coefficient vector does not have the size of the old one");

	cofficients_=newCoef;
}

void combigrid::CombiSchemeBasis::setCoef(int i,double newCoef){
	cofficients_[i]=newCoef;
}
