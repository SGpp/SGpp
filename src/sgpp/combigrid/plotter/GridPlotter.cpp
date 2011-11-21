/*
 * combigrid::GridPlotter.cpp
 *
 *  Created on: Jun 1, 2011
 *      Author: benk
 */

#include "GridPlotter.hpp"

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void combigrid::GridPlotter::plotFullGrid(const std::string& filePath , const FullGridD* fg ,
		std::vector<double>& globalCoord_in , int resolution)
{
	combigrid::Evaluable obj(fg);
	int dim = fg->getDimension();
	plotObject( dim , filePath , &obj , fg->getDomain() , globalCoord_in , resolution);
}

void combigrid::GridPlotter::plotCombiGrid(const std::string& filePath , const combigrid::AbstractCombiGrid* cg ,
		std::vector<double>& globalCoord_in , int resolution)
{
	combigrid::Evaluable obj(cg);
	int dim = cg->getFullGrid(0)->getDimension();
	plotObject( dim , filePath , &obj , cg->getDomain() , globalCoord_in , resolution);
}

void combigrid::GridPlotter::plotObject(int dim ,
		const std::string& filePath ,
		const combigrid::Evaluable* obj ,
		const combigrid::GridDomain* domain ,
		std::vector<double>& globalCoord_in ,
		int resolution)
{
    std::vector<double> result(0);
    std::vector<double> globalCoord = globalCoord_in;
    bool isScaled = false;
    if (domain != 0) {
    	isScaled = domain->get1DDomain(0).isAxisScaled();
    }

	if ( ((domain == 0) || (resolution > 0)) || (!isScaled) )
	{

		if (dim == 1)
		{
		  // if the domain is not existent then just use a regular resolution
		  if (resolution <= 0) resolution = 50;
	      double minX = 0.0;
	      double maxX = 1.0;
	      if (domain!=0) {
	    	  minX = domain->get1DDomain(0).getMinDomain();
	    	  maxX = domain->get1DDomain(0).getMaxDomain();
	      }
	      result.resize( resolution );
	      // loop and evaluate points
	      for (int ii = 0 ; ii < resolution ; ii++){
	    	  globalCoord[0] = minX + (maxX - minX) * (double(ii)/double(resolution-1));
	    	  result[ii] = obj->eval(globalCoord);
	      }

	      // writing file
	      ofstream myfile;
	      myfile.open( filePath.c_str() );
	      myfile << "X = [ " << minX;
	      for (int ii = 1 ; ii < resolution ; ii++){
	    	  myfile << " , " << (minX + (maxX - minX) * (double(ii)/double(resolution-1)));
	      }
	      myfile << "]; \n ";
	      myfile << "res = [ " << result[0];
	      for (int ii = 1 ; ii < resolution ; ii++) {
	    	  myfile << " , " << result[ii];
	      }
	      myfile << "]; \n ";
	      myfile << " plot(X,res); \n " << result[0];
	      myfile.close();
		}
		else
		{
		  double minX = 0.0;
		  double maxX = 1.0;
		  double minY = 0.0;
		  double maxY = 1.0;
		  if (domain!=0) {
		   	  minX = domain->get1DDomain(0).getMinDomain();
		   	  maxX = domain->get1DDomain(0).getMaxDomain();
		   	  minY = domain->get1DDomain(1).getMinDomain();
		   	  maxY = domain->get1DDomain(1).getMaxDomain();
		  }
		  result.resize( resolution * resolution );
		  // loop and evaluate points
		  for (int ii = 0 ; ii < resolution ; ii++){
		   	  for (int jj = 0 ; jj < resolution ; jj++){
		   	     globalCoord[0] = minX + (maxX - minX) * (double(ii)/double(resolution-1));
		   	     globalCoord[1] = minY + (maxY - minY) * (double(jj)/double(resolution-1));
		   	     result[ ii*resolution + jj] = obj->eval(globalCoord);
		   	  }
		   }
		   ofstream myfile;
		   myfile.open(filePath.c_str());
		   myfile << "X = [ " << minX;
		   for (int ii = 1 ; ii < resolution ; ii++){
		       myfile << " , " << (minX + (maxX - minX) * (double(ii)/double(resolution-1)));
		   }
		   myfile << "]; \n ";
		   myfile << "Y = [ " << minY;
		   for (int ii = 1 ; ii < resolution ; ii++){
		       myfile << " , " << (minY + (maxY - minY) * (double(ii)/double(resolution-1)));
		   }
		   myfile << "]; \n ";
		   myfile << "res = [ " << result[0];
		   for (int ii = 1 ; ii < resolution*resolution ; ii++) {
		   	  if ( (ii % resolution ) == 0){
		       	  myfile << " ; " << result[ii];
		      }
		      else{
		   	      myfile << " , " << result[ii];
		      }
		   }
		   myfile << "]; \n ";
		   myfile << "[x,y]=meshgrid(Y,X);\n";
		   myfile << " surf(x,y,res); \n ";
		   myfile.close();
		}
	}
	else
	{
		if (dim == 1)
		{
		      result.resize( domain->get1DDomain(0).axisScaling().size() );
		      // loop and evaluate points
		      for ( unsigned int ii = 0 ; ii < domain->get1DDomain(0).axisScaling().size() ; ii++){
		    	  globalCoord[ 0 ] = domain->get1DDomain(0).axisScaling()[ii];
		    	  result[ii] = obj->eval(globalCoord);
		      }

		       // writing file
		       ofstream myfile;
		       myfile.open(filePath.c_str());
		       myfile << "X = [ " << domain->get1DDomain(0).axisScaling()[0];
		       for ( unsigned int ii = 1 ; ii < domain->get1DDomain(0).axisScaling().size() ; ii++){
		    	   myfile << " , " << domain->get1DDomain(0).axisScaling()[ii];
		       }
		       myfile << "]; \n ";
		       myfile << "res = [ " << result[0];
		       for ( unsigned int ii = 1 ; ii < domain->get1DDomain(0).axisScaling().size() ; ii++) {
		    	   myfile << " , " << result[ii];
		       }
		       myfile << "]; \n ";
		       myfile << " plot(X,res); \n ";
		       myfile.close();
		}
		else
		{
			  result.resize( domain->get1DDomain(0).axisScaling().size() * domain->get1DDomain(1).axisScaling().size() );
			  // loop and evaluate points
			  double res = 0.0;
		      for ( unsigned int ii = 0 ; ii < domain->get1DDomain(0).axisScaling().size() ; ii++){
		    	  for ( unsigned int jj = 0 ; jj < domain->get1DDomain(1).axisScaling().size() ; jj++){
		    	     globalCoord[0] = domain->get1DDomain(0).axisScaling()[ii];
		    	     globalCoord[1] = domain->get1DDomain(1).axisScaling()[jj];
		    	     res = obj->eval(globalCoord);
		    	     result[ ii * domain->get1DDomain(1).axisScaling().size() + jj ] = res;
		    	  }
		      }

		      ofstream myfile;
		      myfile.open(filePath.c_str());
		      myfile << "X = [ " << domain->get1DDomain(0).axisScaling()[0];
		      for ( unsigned int ii = 1 ; ii < domain->get1DDomain(0).axisScaling().size() ; ii++){
		   	     myfile << " , " << domain->get1DDomain(0).axisScaling()[ii];
		      }
		      myfile << "]; \n ";
		      myfile << "Y = [ " << domain->get1DDomain(1).axisScaling()[0];
		      for ( unsigned int ii = 1 ; ii < domain->get1DDomain(1).axisScaling().size() ; ii++){
		   	     myfile << " , " << domain->get1DDomain(1).axisScaling()[ii];
		      }
		      myfile << "]; \n ";
		      myfile << "res = [ " << result[0];
		      for ( unsigned int ii = 1 ; ii < domain->get1DDomain(0).axisScaling().size() *
											   domain->get1DDomain(1).axisScaling().size() ; ii++) {
		    	  if ( (ii % domain->get1DDomain(1).axisScaling().size()) == 0){
		    		  myfile << " ; " << result[ii];
		    	  }
		    	  else{
		   	         myfile << " , " << result[ii];
		    	  }
		      }
		      myfile << "]; \n ";
		      myfile << "[x,y]=meshgrid(Y,X);\n";
		      myfile << " surf(x,y,res); \n ";
		      myfile.close();
		}
	}
}
