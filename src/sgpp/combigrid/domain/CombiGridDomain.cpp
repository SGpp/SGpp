/*
 * CombiGridDomain.cpp
 *
 *  Created on: Apr 4, 2011
 *      Author: benk
 */

#include "CombiGridDomain.hpp"

combigrid::GridDomain::GridDomain(int dim , const std::vector<int>& levels,
        const std::vector<double>& min ,
        const std::vector<double>& max ,
        const combigrid::AbstractStretchingMaker& stretchingMaker) {
	dim_ = dim;
	// add for each dimension
	for (int d = 0; d < dim_ ; d++){
		axisDomains_.push_back( combigrid::Domain1D( levels[d] , min[d] , max[d] , stretchingMaker) );
	}
}

combigrid::GridDomain::GridDomain(int dim , const std::vector< std::vector<double> >& scalings ){
	dim_ = dim;
	// add for each dimension
	for (int d = 0; d < dim_ ; d++){
		axisDomains_.push_back( combigrid::Domain1D( scalings[d] ) );
	}
}

combigrid::GridDomain::GridDomain(int dim ,const std::vector<double>& min ,const std::vector<double>& max ){
	dim_ = dim;
	// add for each dimension
	for (int d = 0; d < dim_ ; d++){
		axisDomains_.push_back( combigrid::Domain1D( min[d] , max[d] ) );
	}
}

void combigrid::GridDomain::transformRealToUnit( std::vector< double >& coords ,
		 const std::vector<int>& levels_in,
		 const std::vector<bool>& boundaryFlag ) const {
	// for each dimension make the transformation
	//int verb = 6;
	double tmp = 0.0;
	//COMBIGRID_OUT_LEVEL3( verb , " combigrid::GridDomain::transformRealToUnit() ");
	for (int d = 0; d < dim_ ; d++){
		axisDomains_[d].transformRealToUnit( coords[d] , tmp , levels_in[d] , boundaryFlag[d] );
		coords[d] = tmp;
	}
}
void combigrid::GridDomain::printDomain(){
	std::cout<<"------------------"<<std::endl;
	for (int d=0;d<dim_;d++){
		std::cout<<axisDomains_[d].getMinDomain()<<"\t"<<axisDomains_[d].getMaxDomain()<<"\t"<<axisDomains_[d].isAxisScaled()<<std::endl;
	}
	std::cout<<"------------------"<<std::endl;
}


