/*
 * CombiS_CT.hpp
 *
 *  Created on: Feb 26, 2011
 *      Author: benk
 */

#ifndef COMBIS_CT_HPP_
#define COMBIS_CT_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "CombiCombiSchemeBasis.hpp"

namespace combigrid {

/** class to model the standard combination technique. <br>
 * For grids with boundary there will be only the trapezoid grid.
 * The linear grid with boundary has too many boundary points
 * and is not often used in practical applications.
 *  */
class S_CT : public CombiSchemeBasis{

public:


	/** simple Ctor with the same level */
	S_CT( int dim , int level ) : CombiSchemeBasis(dim) {

		std::vector< int > levels( dim_ , level);
		// the rations for the dimension adaptive case
		std::vector< double > ratio_( dim_ , 1.0 );
		// the truncation level, if there is none then l_user = 1 for all dimensions ( trapezoidal ) for l_user = 0 (linear)
		std::vector< int > l_user_( dim_ , 1 );

		// call the init function
		init( levels , ratio_ , l_user_ );
	}


	/** simple Ctor with the same level and level of truncation*/
	S_CT( int dim , int level , int level_truncation) : CombiSchemeBasis(dim) {

		std::vector< int > levels( dim_ , level);
		// the rations for the dimension adaptive case
		std::vector< double > ratio_( dim_ ,1.0);
		// the truncation level, if there is none then l_user = 1 for all dimensions ( trapezoidal ) for l_user = 0 (linear)
		std::vector< int > l_user_( dim_ , level_truncation );

		// call the init function
		init( levels , ratio_ , l_user_ );
	}


	/** S-CT Ctor with dimension adaptive levels */
	S_CT( int dim , const std::vector< int >& levels ) : CombiSchemeBasis(dim) {

		std::vector< int > levels_tmp = levels;
		// the rations for the dimension adaptive case
		std::vector< double > ratio_( dim_ ,1.0);
		// the truncation level, if there is none then l_user = 1 for all dimensions ( trapezoidal ) for l_user = 0 (linear)
		std::vector< int > l_user_( dim_ , 1 );

		// call the init function
		init( levels_tmp , ratio_ , l_user_ );
	}


	/** S-CT Ctor with dimension homogenious levels but with dimension adaptive truncation */
	S_CT( int dim , int level , const std::vector< int >& levels_trunk ) : CombiSchemeBasis(dim) {

		std::vector< int > levels( dim_ , level );
		// the rations for the dimension adaptive case
		std::vector< double > ratio_( dim_ ,1.0);
		// the truncation level, if there is none then l_user = 1 for all dimensions ( trapezoidal ) for l_user = 0 (linear)
		std::vector< int > l_user_;

		// this is trapezoid grid
		l_user_ = levels_trunk;

		// call the init function
		init( levels , ratio_ , l_user_ );
	}


	/** S-CT Ctor with dimension adaptive levels and with dimension adaptive levels and truncation */
	S_CT( int dim , const std::vector< int >& levels , const std::vector< int >& levels_trunk ) : CombiSchemeBasis(dim) {

		std::vector< int > levels_tmp = levels;
		// the rations for the dimension adaptive case
		std::vector< double > ratio_( dim_ ,1.0);
		// the truncation level, if there is none then l_user = 1 for all dimensions ( trapezoidal ) for l_user = 0 (linear)
		std::vector< int > l_user_ = levels_trunk;

		// call the init function
		init( levels_tmp , ratio_ , l_user_ );
	}


protected:


	/** the init function which is general for all types of combi schemes
	 * */
	void init( std::vector< int >& levels , std::vector< double >& ratio_ ,	std::vector< int >& l_user_ ) {

		std::vector< int > v(0);

		// the ratio for dimension adaptivity
		ratio_.resize( dim_ , 1.0 );

		// get the maximum level and get the
		int n = levels[0] , max = 0;
		for (int i = 1 ; i < dim_ ; i++ ){
			if ( levels[i] > n )
				n = levels[i];
			if ( l_user_[i] > max )
				max = l_user_[0];
		}
		for (int i = 0 ; i < dim_ ; i++ ){
			ratio_[i] = (double)n/(double)levels[i];
		}
		n = n - max;

		// add the different full grids
		double combi = 0.0;
		for (int d = 0 ; d < dim_ ; d++ )
		{
				combi = combigrid::combination( dim_-1 , d );
				int oldsize = levels_vector_.size();
				if ( d%2 != 0 ) combi=-combi;
				// call the recursive function to add the spaces
				getTrapezoidsums( v , dim_ , n-d , ratio_ , l_user_);
				// add the coefficients
				for (int i = oldsize ; i < (int)levels_vector_.size() ; i++ )
					cofficients_.push_back(combi);
		 }
		 // remove the duplicated spaces
		 removeDuplicates();
	}


	/** recursive function for the combination scheme
	 * */
	void getTrapezoidsums(std::vector<int>& v , size_t dim , int sum , std::vector< double >& ratio_ ,	std::vector< int >& l_user_)
	{
		/* Takes recursively every possible combination of numbers which add up to sum creating a linear boundary grid for each one
		 * The levels of the full grids must be greater than l_user*/
		// code from Aliz
		if (dim==1)
		{
			v.push_back( int(sum/ratio_[v.size()]) + l_user_[v.size()] );
			//add v to the level vectors
			levels_vector_.push_back(v);
			v.pop_back();
		}
		else{
			for (int i=0 ; i <= sum ; i++ )
			{
				v.push_back( (int)(i/ratio_[v.size()]) + l_user_[v.size()] );
				getTrapezoidsums( v , dim-1 , sum-i , ratio_ , l_user_ );
				v.pop_back();
			}
		}
	}

};
}


#endif /* COMBIS_CT_HPP_ */
