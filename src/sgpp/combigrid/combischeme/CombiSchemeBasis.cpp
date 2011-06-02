/*
 * CombiCombiSchemeBasis.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */


#include "combigrid/combischeme/CombiSchemeBasis.hpp"

using namespace std;


void combigrid::CombiSchemeBasis::removeDuplicates(){
	// this is Alize's code 1 to 1
	size_t i=0;
	size_t j;
	std::vector< std::vector<int> >::iterator gt = levels_vector_.begin();
	std::vector<double>::iterator ct = cofficients_.begin();
	int size_fg = levels_vector_.size();
	while( (int)i < size_fg-1 ){
		j = i + 1;
		while( (int)j < size_fg )
		{
			if ( levels_vector_[i] == levels_vector_[j] ){
				levels_vector_.erase(gt+j);
				cofficients_[i]+=cofficients_[j];
				cofficients_.erase(ct+j);
				size_fg--;
			}
			else j++;
		}
		if ( cofficients_[i] == 0.0 ){
			levels_vector_.erase(gt+i);
			cofficients_.erase(ct+i);
			size_fg--;
			if ( i == 0 ){
				gt = levels_vector_.begin();
				ct = cofficients_.begin();
			}
		}
		else i++;
	}
}
