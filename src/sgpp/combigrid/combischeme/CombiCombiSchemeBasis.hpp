/*
 * CombiCombiSchemeBasis.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef COMBICOMBISCHEMEBASIS_HPP_
#define COMBICOMBISCHEMEBASIS_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"

using namespace std;

namespace combigrid {

	/** base class for any combi scheme. From this class should all the scheme classes derived e.g. S-CT, TS-CT, ... <br>
	 * */
	class CombiSchemeBasis {

	public:

		/** Ctor with no argument */
		CombiSchemeBasis(int dim) { dim_ = dim; levels_vector_.resize(0); cofficients_.resize(0); }

		/** return the dimension of the combi scheme */
		inline int getDim() const { return dim_; }

		/** number of subsapces */
		inline int getNrSapces() const { return levels_vector_.size(); }

		/** returns the level vector for one subspace */
		inline const std::vector<int>& getLevel(int i) const { return levels_vector_[i]; }

		/** returns the coefficient for one subspace */
		inline double getCoef(int i) const { return cofficients_[i]; }

	protected:

		/** function which removes the duplicate spaces */
		void removeDuplicates(){
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

		/** the dimension of the scheme */
		int dim_;

		/** the level vector for each space */
		std::vector< std::vector<int> > levels_vector_;

		/** the coefficients for the spaces */
		std::vector< double > cofficients_;
	};

}

#endif /* COMBICOMBISCHEMEBASIS_HPP_ */
