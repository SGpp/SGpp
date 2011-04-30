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

		/** Ctor */
		CombiSchemeBasis(int dim, int level) {
			dim_ = dim; levels_vector_.resize(0); levels_.resize(dim_,level); cofficients_.resize(0);
		}

		/** Ctor */
		CombiSchemeBasis(int dim, const std::vector<int>& levels) {
			dim_ = dim; levels_vector_.resize(0); levels_ = levels; cofficients_.resize(0);
		}

		/** return the dimension of the combi scheme */
		inline int getDim() const { return dim_; }

		/** number of subsapces */
		inline int getNrSapces() const { return levels_vector_.size(); }

		/** returns the level vector for one subspace */
		inline const std::vector<int>& getLevel(int i) const { return levels_vector_[i]; }

		/** the levels of the full grid which we extrapolate with the combi scheme*/
		inline const std::vector<int>& getMaxLevel() const { return levels_; }

		/** returns the coefficient for one subspace */
		inline double getCoef(int i) const { return cofficients_[i]; }

	protected:

		/** function which removes the duplicate spaces */
		void removeDuplicates();

		/** the dimension of the scheme */
		int dim_;

		/** the level vector for each space */
		std::vector< std::vector<int> > levels_vector_;

		/** the levels of the full grid which we*/
		std::vector<int> levels_;

		/** the coefficients for the spaces */
		std::vector< double > cofficients_;
	};

}

#endif /* COMBICOMBISCHEMEBASIS_HPP_ */
