/*
 * combigrid::RunTikhonov.cpp
 *
 *  Created on: May 30, 2011
 *      Author: benk
 */

#include "RunTikhonov.hpp"
#include "solver/multigridFG/operators/TikhonovOperator.hpp"
#include "solver/multigridFG/multigrid/Multigrid.hpp"


void combigrid::RunTikhonov::readInInput(
		const std::string& XfileS ,
		const std::string& YFileS ,
		int& dimensions ,  int& nrPoints ,
		std::vector<double>& XCoords,
		std::vector<double>& YPoint) {

	int verb = 3;
    dimensions = 0; nrPoints = 0;

    std::ifstream Xfile;
    std::ifstream Yfile;

    Xfile.open( XfileS.c_str() , std::ios::in );
    Yfile.open( YFileS.c_str() , std::ios::in );

    // the string which we read in one line
    std::string str_tmp;
    double d1;
    int i;
    char *pEnd,*pStart;

    // read in one line to detect the dimension
    getline( Xfile , str_tmp );
    dimensions = 0;
    pStart = &(str_tmp[0]); //str_tmp.c_str();
    std::cout << " first line: " << str_tmp << std::endl;
    while ( pStart[0] != 0 ){
        d1 = strtod (pStart , &pEnd);
        std::cout << " number: " << d1 << " pEnd:" << (int)(pEnd[0]) << std::endl;
        pStart = pEnd;
        if (pStart != NULL){
            dimensions++;
            XCoords.push_back(d1);
        }
    }
    std::cout << " dim: " << dimensions << std::endl;

    // now we have the dimension, read in the X
    getline( Xfile , str_tmp );
    while (!Xfile.eof())
    {
       // for each line get the
       pStart = &(str_tmp[0]); //str_tmp.c_str();
       for(i = 0 ; i < dimensions ; i++){
          d1 = strtod (pStart , &pEnd);
          pStart = pEnd;
          XCoords.push_back(d1);
       }
       getline( Xfile , str_tmp );
    }

    // read in Y
    getline( Yfile , str_tmp );
    while (!Yfile.eof())
    {
       // there is only one number per line
       pStart = &(str_tmp[0]); //str_tmp.c_str();
       d1 = strtod (pStart , &pEnd);
       YPoint.push_back(d1);
       getline( Yfile , str_tmp );
       nrPoints++;
    }

    // print the results
    int j = 0;
    if (verb > 3){
    	for( i = 0 ; i < nrPoints ; i++){
    		std::cout << " X["<<i<<"] ( " << XCoords[i*dimensions];
    		for ( j = 1 ; j < dimensions ; j++) {
    			std::cout << "," << XCoords[i*dimensions+j];
    		}
    		std::cout << " ) Y = " << YPoint[i] << std::endl;
    	}
    }
}

FullGridD* combigrid::RunTikhonov::computeFGTikhonov(
		combigrid::GridDomain& domain ,
		const std::vector<int>& levels,
		double lambda,
		const std::string& Xfile ,
		const std::string& YFile ){

	int dimensions = 0;
	int nrPoints = 0;
	int verb = 6;
	std::vector<double> XCoords;
	std::vector<double> YPoint;
	std::vector<double> unknowns;
	FullGridD* fg;
	// first read in the results
	COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... read INPUT");
	combigrid::RunTikhonov::readInInput(Xfile , YFile , dimensions , nrPoints , XCoords, YPoint);

	// create the full grid
	COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create Full Grid");
	fg = new FullGridD( dimensions , levels );
	fg->createFullGrid();
	fg->setDomain( &domain );

	// create the Tikhonov operator
	COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create combigrid::TikhonovOperator");
	combigrid::TikhonovOperator op( fg , nrPoints, lambda , &XCoords , &YPoint );

	// create the multigrid method
	COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create combigrid::Multigrid");
	combigrid::Multigrid multigrid( &(op) , fg );

	// solve using only the smoother
	COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... solve multigird");
	unknowns.resize(fg->getNrElements(),0.0);
	multigrid.solveCS( unknowns , 1e-8 , false );
	//multigrid.solveSmoothing( unknowns , 1e-8);
	//multigrid.solveCG( unknowns , 1e-6 );
	//multigridFAS.solveFAS( unknowns , 1e-8 );

	// copy the solution back
	COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... write solution back and return");
	for (int i = 0 ; i < fg->getNrElements() ; i++){
		fg->getElementVector()[i] = unknowns[i];
		//COMBIGRID_OUT_LEVEL2( verb , "unknowns["<<i<<"]="<<unknowns[i]);
	}

	// return full grid
	return fg;
}
