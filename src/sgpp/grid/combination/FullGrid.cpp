#include "FullGrid.hpp"
#include "LinearFullGrid.hpp"
using namespace std;
namespace sg{


    FullGrid* FullGrid::createLinearFullGrid(size_t dim, vector<level_t> *inlevel)
    {
        return new LinearFullGrid(dim,inlevel);
    }


    FullGrid* FullGrid::createLinearBoundaryFullGrid(size_t dim, vector<level_t> *inlevel)
    {
        return new FullGrid(dim,inlevel);
    }


   void FullGrid::getCoords(size_t index,DataVector &v)
    {
        int ind=0;
        int aux=index;
        /**
         *The number of gridpoints on every direction is Ni:=2^level[i]+1, so the index of a gridpoint in a given direction
         *can be computed by the formula ind(i)=[index/(N0*N1...        *Ni-1)] mod Ni, where index is the position of the gridpoint
         *in the vector.
         *The coordinate of a gridpoint in direction 'i' is ind(i)*2^(-level[i]), considering the grid of dimension 1.0
         */
        for (size_t i=0;i<dim;i++){
               ind=aux%((int)powOfTwo[level[i]]+1);
	           aux=aux/((int)powOfTwo[level[i]]+1);
               if (ind==0) v[i]=0;
               else
               v[i]=((double)ind)*negPowOfTwo[level[i]];
        }
        if (FullGrid::boundingBox!=0)
        {
        	 for (size_t i=0;i<dim;i++)
        		 v[i]=v[i]*FullGrid::boundingBox->getIntervalWidth(i)+FullGrid::boundingBox->getIntervalOffset(i);
        }
    }


   std::string FullGrid::getCoordsString(size_t index)
    {
    	std::stringstream return_stream;

    	// switch on scientific notation:
		//return_stream << std::scientific;
        DataVector coords(dim);
        getCoords(index,coords);
        //return_stream.precision(15);
    	for(size_t i = 0; i < dim; i++)
    	{
    		return_stream << std:: scientific<<coords[i];    		
    		if (i < dim-1)
    		{
    			return_stream << " ";
    		}
    	}

    	return return_stream.str();
    }


   double FullGrid::eval(DataVector& p)
   {
  	 int ii,i,tmp_val,vv,nr;
  	 int jj;
  	 double dist,baseVal,ofs,len;
  	 double normcoord;
  	 /**If the boundingBox of the fullgrid is null we consider the default [0,1]^d hypercube*/
  	 if (boundingBox!=0){
  		 for ( ii = dim-1 ; ii >=0; ii--){
  			 ofs=boundingBox->getIntervalOffset(ii);
				 len=boundingBox->getIntervalWidth(ii);
  		     // scale to the reference [0,1] interval the intersection point
  		     normcoord=(p[ii]-ofs)*powOfTwo[level[ii]]/len;
  		     aindex[ii] = floor(normcoord);
  		     aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
			 aindex[ii] = (aindex[ii] >= (int)(this->length(ii) - 1)) ? (int)(this->length(ii)-2) : aindex[ii];
  		     dist = aindex[ii]+1-normcoord;
  		     intersect[2*ii] = dist;
  		     intersect[2*ii+1] = 1- dist;
  		  }
  	 }
  	 else
  	 {
			 for ( ii = dim-1 ; ii >=0; ii--){
			     // scale to the reference [0,1] interval the intersection point
				 normcoord=p[ii]*powOfTwo[level[ii]];
				 aindex[ii]=floor(normcoord);
				 aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
				 aindex[ii] = (aindex[ii] >= (int)(this->length(ii) - 1)) ? (int)(this->length(ii)-2) : aindex[ii];
				 dist = aindex[ii]+1-normcoord;
				 intersect[2*ii] = dist;
				 intersect[2*ii+1] = 1 - dist;
			 }
  	 }
  	 value = 0.0;
  	 nr=powOfTwo[dim];
  	 for (ii=0; ii < nr; ii++){
  	        baseVal = 1;//we will store here the coefficient for the next corner point of the cuboid
  	        i = 0;
  	        tmp_val = ii;
  	        vv = 0;
  	        // compute the "dim" dimensional basis function value(for one node) and the corresponding vector index
  	        for (jj = dim-1 ; jj >=0; jj--){
  	              vv = (tmp_val & 1); // tmp % 2 ;
  	              baseVal = baseVal * intersect[2*jj+vv];
  	              i = i*this->length(jj)+aindex[jj]+vv;
  	              tmp_val = tmp_val >> 1; //tmp_val / 2;
  	              //std::cout << "FullGrid::eval jj:" << jj << " , this->length(jj):" << this->length(jj) << " , i:" << i
  	            //		  << ", aindex[jj]:" << aindex[jj] << " , vv:" << vv<< std::endl;
  	        }
  	        // multiply the basis function value with the coefficient
  	        //if (baseVal!=0) //todo:
  	        value += baseVal * vec[i];
  	  }
  	  // just return the value
  	  return value;
   }
}
