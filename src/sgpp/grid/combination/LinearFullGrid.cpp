/*
 * LinearFullGrid.cpp
 *
 *  Created on: Jan 8, 2011
 *      Author: benk
 */



#include "LinearFullGrid.hpp"

#include <string.h>
#include <stdlib.h>
using namespace sg::base;

namespace sg{

double LinearFullGrid::eval(DataVector& p)
{
	    int ii,i,tmp_val,vv,nr;
	    int jj;
	    double dist,baseVal,ofs,len;
	    double normcoord;
	    double ret_val = 0.0;
		double intersect[2*dim];
		int aindex[dim];
	    /**If the boundingBox of the fullgrid is null we consider the default [0,1]^d hypercube*/
	    if (boundingBox!=0){
	    		 for ( ii = dim-1 ; ii >=0; ii--){
	    			 ofs=boundingBox->getIntervalOffset(ii);
					 len=boundingBox->getIntervalWidth(ii);
	    		     // scale to the reference [0,1] interval the intersection point
	    		     normcoord=(p[ii]-ofs)*powOfTwo[level[ii]]/len;
	    		     aindex[ii]=floor(normcoord-1e-14);
	    		     dist = aindex[ii]+1-normcoord;
					 intersect[2*ii] = dist;
					 intersect[2*ii+1] = 1 - dist;
	    		     aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
	    		     aindex[ii] = (aindex[ii] >= (int)(this->length(ii))) ? (int)(this->length(ii)) : aindex[ii];
	    		     // for special case make some adjustment
		             if (aindex[ii]==0){
		            	  aindex[ii] = 1;
		            	  intersect[2*ii] = (intersect[2*ii]+1); //extrapolated basis value
		            	  intersect[2*ii+1] = (intersect[2*ii+1]-1); //extrapolated basis value
		              }
		              else{
		               // make extrapolation at the last cell (where there is no boundary point)
		               if (aindex[ii] >= (int)this->length(ii)){
		                  aindex[ii] = (int)this->length(ii)-1;
		            	  intersect[2*ii] = (intersect[2*ii]-1); //extrapolated basis value
		            	  intersect[2*ii+1] = (intersect[2*ii+1]+1);  //extrapolated basis value
		               }
		             }
	    		  }
	    }
	    else
	    {
				 for ( ii = dim-1 ; ii >=0; ii--){
				     // scale to the reference [0,1] interval the intersection point
					 normcoord=p[ii]*powOfTwo[level[ii]];
					 aindex[ii]=floor(normcoord-1e-14);
					 dist = aindex[ii]+1-normcoord;
					 intersect[2*ii] = dist;
					 intersect[2*ii+1] = 1 - dist;
	    		     aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
	    		     aindex[ii] = (aindex[ii] >= (int)(this->length(ii))) ? (int)(this->length(ii)) : aindex[ii];
	    		     // for special case make some adjustment
		             if (aindex[ii]==0){
		            	  aindex[ii] = 1;
		            	  intersect[2*ii] = (intersect[2*ii]+1); //extrapolated basis value
		            	  intersect[2*ii+1] = (intersect[2*ii+1]-1); //extrapolated basis value
		              }
		              else{
		               // make extrapolation at the last cell (where there is no boundary point)
		               if (aindex[ii] >= (int)this->length(ii)){
		                  aindex[ii] = (int)this->length(ii)-1;
		            	  intersect[2*ii] = (intersect[2*ii]-1); //extrapolated basis value
		            	  intersect[2*ii+1] = (intersect[2*ii+1]+1);  //extrapolated basis value
		               }
		             }
				 }
	    }

	    size_t d=dim;
	    // detect the "real" dimensions, where there is only one point in one dimension
	    for (ii=0;ii < (int)dim;ii++){
		  if (level[ii]==1)
		  {
			 d--;
	      }
	    }

	    nr=powOfTwo[d];
	    // calculate the
	    for (ii=0; ii < nr; ii++){
	        baseVal = 1.0;
	        i = 0;
	        tmp_val = ii;
	        vv = 0;
	        for (jj = dim-1 ; jj >=0; jj--){
	           // skip the dimensions where there in only one point
	           if (this->level[jj]!=1)
	           {
	              vv = (tmp_val & 1);
	              // make interpolation (as in the usual case)
	              i = i*this->length(jj)+aindex[jj]-1+vv;
	              baseVal = baseVal * intersect[2*jj+vv];
	              //
	              tmp_val = tmp_val >> 1; //tmp_val / 2;
	           }
	        }// end for loop

	        // multiply the basis function value with the coefficient
	        ret_val += baseVal* vec[i];

	    }
	    value = ret_val;
	    // just return the value
	    return ret_val;
}

}
