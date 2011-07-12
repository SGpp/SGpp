/*
 * CombiDomain1D.cpp
 *
 *  Created on: Apr 4, 2011
 *      Author: benk
 */

#include "CombiDomain1D.hpp"
#include <vector>
#include <algorithm>



combigrid::Domain1D::Domain1D(double min, double max) {
	isStretched_ = false;
	level_ = -1;
	min_ = min;
	max_ = max;
}


combigrid::Domain1D::Domain1D(const std::vector<double>& inputStretching){
	isStretched_ = true;
	level_ = ::round( ::log((double)inputStretching.size())/ ::log(2.0) );

	// test if the vector size match
	COMBIGRID_ERROR_TEST( (int)inputStretching.size() == combigrid::powerOfTwo[level_]+1 ,
			"Domain1D::Domain1D , Input vector size must be 2^L+1 , inputStretching.size():" << inputStretching.size()
			<< " , level_:" << level_);

	// copy the vector and sort
	stretching_ = inputStretching;
	std::sort( stretching_.begin(),stretching_.end() );
	min_ = stretching_[0];
	max_ = stretching_[stretching_.size()-1];
}


combigrid::Domain1D::Domain1D(int level, double min, double max, const combigrid::AbstractStretchingMaker& stretching){
	isStretched_ = true;
	level_ = level;
	stretching.get1DStretching( level , min , max , stretching_);
	min_ = stretching_[0];
	max_ = stretching_[stretching_.size()-1];
}


void combigrid::Domain1D::transformRealToUnit(double coordReal,
		 double& coordUnit ,
		 int level_in ,
		 bool noBoundary ) const{

	//int verb = 6;

	if (isStretched_){
		int startInd = 0 ,mid = 0;
		int endInd = stretching_.size() - 1;
		double intersec = 0.0;
		int level_diff = (level_ < level_in)? 0 : level_ - level_in;
		// stop when the difference is one, which means we found the cell
		while ( endInd - startInd > combigrid::powerOfTwo[level_diff] ) {
			mid = ((endInd + startInd)/2);
			// make the bisection
			if ( stretching_[mid] < coordReal){
				startInd = mid;
			} else {
				endInd = mid;
			}
		}

		// startInd should be now at the beginning of the cell
		intersec = (stretching_[endInd] - coordReal) /
				(stretching_[endInd] - stretching_[startInd]);
		//   - make the transformation to the unit domain -> unitCoords
		coordUnit = ( ((double)endInd)*combigrid::oneOverPowOfTwo[level_diff] - intersec) /
				combigrid::powerOfTwo[ level_-level_diff];

		if ( noBoundary ){
			int offs = combigrid::powerOfTwo[level_diff];
			// for boundary cells make extrapolation
			if (startInd == 0) {
				double h = stretching_[endInd+offs] - stretching_[startInd+offs];
				intersec = (stretching_[endInd] - coordReal)/h;
				coordUnit = ( ((double)endInd)*combigrid::oneOverPowOfTwo[level_diff] - intersec) /
						combigrid::powerOfTwo[ level_-level_diff];
				//COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::tran 1");
			}
			if (endInd == (int)stretching_.size() - 1){
				double h = stretching_[startInd] - stretching_[startInd-offs];
				intersec = (coordReal - stretching_[startInd])/h;
				coordUnit = ( ((double)(endInd-offs))*combigrid::oneOverPowOfTwo[level_diff] + intersec) /
						combigrid::powerOfTwo[ level_-level_diff];
				//COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::tran 2 , intersec:" << intersec << " , endInd:"<<endInd);
				//COMBIGRID_OUT_LEVEL3( verb , " stretching_[startInd]:" << stretching_[startInd]<< " , stretching_[endInd]:" << stretching_[endInd]);
				//COMBIGRID_OUT_LEVEL3( verb , " h:" << h << " , h_old:" << stretching_[endInd] - stretching_[startInd] );
			}
		}
		//COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::transformRealToUnit coordReal:" << coordReal << " coordUnit:"
		//		<< coordUnit << "  level_in:"<<level_in);
	}
	else
	{ // no stretching , just simple scaling
		coordUnit = (coordReal - min_)/(max_ - min_);
		//COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::transformRealToUnit NO STRETCH coordReal:" << coordReal << " coordUnit:"
		//		<< coordUnit << "  level_in:"<<level_in);
	}
}


void combigrid::Domain1D::transformUnitToReal( int level , int index ,
		double& realCoord) const {
	if (isStretched_){
		// get the stretched index
		int level_diff = (level_ < level) ? 0 : level_ - level;
		realCoord = stretching_[ index * combigrid::powerOfTwo[level_diff] ];
	}
	else
	{   // no stretching , just simple scaling
		realCoord = min_ + (max_ - min_)*((double)index)*oneOverPowOfTwo[level];
	}
}


void combigrid::Domain1D::findEntry(double coordReal, int level_in ,
		int& startIndex , double& intersect) const {
	if (isStretched_)
	{
		int startInd = 0 ,mid = 0;
		int endInd = stretching_.size() - 1;
		int level_diff = (level_ < level_in)? 0 : level_ - level_in;
		// stop when the difference is one, which means we found the cell
		while ( endInd - startInd > combigrid::powerOfTwo[level_diff] )
		{
			mid = ((endInd + startInd)/2);
			// make the bisection
			if ( stretching_[mid] < coordReal){
				startInd = mid;
			} else {
				endInd = mid;
			}
		}
		// startInd should be now at the beginning of the cell
		intersect = ( coordReal - stretching_[startInd] ) / //stretching_[endInd] -
				(stretching_[endInd] - stretching_[startInd]);
		// this must be the start index of the grid, not from the stretching (the levels could be different)
		startIndex = startInd/combigrid::powerOfTwo[level_diff];
	}
	else
	{
		// for the non-stretched case
		double unitC = (coordReal - min_)/(max_ - min_);
		startIndex = ::floor( (double)combigrid::powerOfTwo[level_in]*unitC );
		startIndex = (startIndex < 0) ? 0 : startIndex;
		startIndex = (startIndex >= combigrid::powerOfTwo[level_in]-1) ? combigrid::powerOfTwo[level_in]-1 : startIndex;
		intersect =  (unitC*(double)combigrid::powerOfTwo[level_in]  - (double)(startIndex));
	}
}


void combigrid::Domain1D::getMeshWidth(int index , int level_in , double& h0 , double& h1) const {
	if (isStretched_)
	{
		int level_diff = (level_ < level_in)? 0 : level_ - level_in;
		// index checking should be done in the debug mode
		h0 = stretching_[ index*combigrid::powerOfTwo[level_diff] ] - stretching_[ (index-1)*combigrid::powerOfTwo[level_diff] ];
		h1 = stretching_[ (index+1)*combigrid::powerOfTwo[level_diff] ] - stretching_[ index*combigrid::powerOfTwo[level_diff] ];
	}
	else
	{
		h1 = h0 = (max_ - min_) * ( 1 / (double)combigrid::powerOfTwo[level_in] );
	}
}

