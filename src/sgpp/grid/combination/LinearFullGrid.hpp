#ifndef LINEARFULLGRID_HPP
#define LINEARFULLGRID_HPP

#include "FullGrid.hpp"

namespace sg{

/**
  *A  linear fullgrid without boundary
  *Every dimension has 2^level[i]-1 gridpoints
  */
class LinearFullGrid: public FullGrid {
public:
/**
 *Constructor for the linear fullgrid without boundary with 2^level[i]-1 points in every direction
 *@param indim the dimension of the fullgrid
 *@param inlevel the level vector for the fullgrid
*/
LinearFullGrid(size_t indim, vector<level_t> *inlevel): FullGrid(indim){
    for (size_t i=0;i<FullGrid::dim;i++){
		FullGrid::level[i]=inlevel->at(i);
		FullGrid::size=FullGrid::size*(powOfTwo[inlevel->at(i)]-1);
    }
    vec.resize(size,0);
    intersect=new double[2*dim];
    aindex=new int[dim];
}
/**
 *Returns the gridpoint for given index vector 
 *@param index the indexes of the gridpoint in every direction
 *@return acces to the gridpoint value in the data vector
 */
double& at(index_t *index){
    size_t ind=0;
    for (int i=dim-1;i>=0;i--)
    	    	ind=(powOfTwo[level[i]]-1)*ind+index[i];
    return vec.at(ind);
}
/**
 *Returns the coordinate vector for a given gridpoint
 *@param index the index of the gridpoint in the data vector
 *@param v the datavector in which the coordinates will be placed
 */
void getCoords(size_t index,DataVector &v)
    {
        int ind=0;
        int aux=index;
        /**
         *The number of gridpoints on every direction is Ni:=2^level[i]-1, so the index of a gridpoint in a given direction
         *can be computed by the formula ind(i)=[index/(N0*N1...        *Ni-1)] mod Ni, where index is the position of the gridpoint
         *in the vector.
         *The coordinate of a gridpoint in direction 'i' is (ind(i)+1)*2^(-level[i]), considering the grid of dimension 1.0
         */
        for (size_t i=0;i<dim;i++){              
               ind=(aux%(powOfTwo[level[i]]-1));
	       aux=aux/(powOfTwo[level[i]]-1);
               v[i]=((double)(ind+1))*negPowOfTwo[level[i]];
        }
        if (FullGrid::boundingBox!=0)
        {
        	 for (size_t i=0;i<dim;i++)
        		 v[i]=v[i]*FullGrid::boundingBox->getIntervalWidth(i)+FullGrid::boundingBox->getIntervalOffset(i);
        }
     }
inline double getCoord(size_t d,size_t index){
    int ind=index;
    double c=0;
    /**
    *Computes coordinate of a gridpoint in direction d as: [[index/(N0*N1*N2...*Nd-1)] %Nd+1] *2^(-level[d])
    *Ni=2^level[i]+1
    */
    for (size_t i=0;i<d;i++)
    	   {
    		 ind=ind / (powOfTwo[level[i]]-1);
    	   }
    	   ind=ind%(powOfTwo[level[d]]-1);
           c=((double)(ind+1))*negPowOfTwo[level[d]];
    if (FullGrid::boundingBox!=0)
          	        {
          	        		 c=c*FullGrid::boundingBox->getIntervalWidth(d)+FullGrid::boundingBox->getIntervalOffset(d);
          	        }
    return c;
}
/**
 * Returns the number of gridpoints in direction d which is 2^level[d]-1
 */
std::string getCoordsString(size_t index)
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
size_t length(size_t d)
{
   return powOfTwo[FullGrid::level[d]]-1;
}
/**
 *Returns the starting index of the gridpoints 
 */
inline size_t startindex()
     {
        return 1;
      }
/**
 * Returns the type of the fullgrid
 * */
const char* getType()
{
      return "linear";
}
/**
 * Evaluates the value of the function in a given point
 * */
double eval(DataVector& p)
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
   	    		     		 aindex[ii]=floor(normcoord);
   	    		     	     dist= aindex[ii]+1-normcoord;
   	    		     	     intersect[2*ii] =dist;
   	    		     	     intersect[2*ii+1] =1- dist;
   	    		     	    }
   	    	 }
   	    	 else
   				 for ( ii = dim-1 ; ii >=0; ii--){
   				// scale to the reference [0,1] interval the intersection point
   					 normcoord=p[ii]*powOfTwo[level[ii]];
   					 aindex[ii]=floor(normcoord);
   					 dist= aindex[ii]+1-normcoord;
   					 intersect[2*ii] =dist;
   					 intersect[2*ii+1] =1- dist;
   					}

   	 value = 0.0;
   	 size_t d=dim;
   	 for (ii=0;ii<dim;ii++)
   		 if (level[ii]==1)
   			 d--;
   	 nr=1<<d;
   	 for (ii=0; ii < nr; ii++){
   	        baseVal = 1.0;
   	        i = 0;
   	        tmp_val = ii;
   	        vv = 0;
   	        for (jj = dim-1 ; jj >=0; jj--){   	              // tmp % 2 ;
//   	              if ((this->level[jj]==1))
//   	            	  {
//						  i=i*this->length(jj);
//						  baseVal = baseVal * intersect[2*jj+vv];
//   	            	  }
//   	              else
//   	              if ((aindex[jj]==0)&&(vv==0)){
//   	            	  i = i*this->length(jj)+aindex[jj];
//   	            	  baseVal*=2* intersect[2*jj+vv];
//   	              }
//   	              else
//   	               if ((aindex[jj]==this->length(jj))&&(vv==1)){
//   	            	 i = i*this->length(jj)+aindex[jj]-1;
//   	            	 baseVal*=2* intersect[2*jj+vv];
//   	               }
//   	              else{
//   	            	 i = i*this->length(jj)+aindex[jj]-1+vv;
//   	            	 baseVal = baseVal * intersect[2*jj+vv];
//   	              }
   	              if (this->level[jj]!=1)
   	            	  {
   	            	vv = (tmp_val & 1);
   	              if (aindex[jj]==0){
   	            	  i = i*this->length(jj)+aindex[jj]+vv;
   	            	  if (vv==0)
   	            		  baseVal*= (intersect[2*jj+vv]+1);
   	            	  else
   	            		  baseVal*=(intersect[2*jj+vv]-1);
   	              }
   	              else
   	               if (aindex[jj]==this->length(jj)){
   	            	 i = i*this->length(jj)+aindex[jj]-1;
   	            	 if (vv==0)
   	            	   	  baseVal*= (intersect[2*jj+vv]-1);
   	            	 else
   	            	   	  baseVal*=(intersect[2*jj+vv]+1);
   	               }
   	              else{
   	            	 i = i*this->length(jj)+aindex[jj]-1+vv;
   	            	 baseVal = baseVal * intersect[2*jj+vv];
   	              }

   	              tmp_val = tmp_val >> 1; //tmp_val / 2;
   	        }}
   	        // multiply the basis function value with the coefficient
   	        if (baseVal!=0)
   	        	value += baseVal* vec[i];

   	    }
   	        // just return the value
   	    return value;
       return 0;
    }
};
}
#endif /* LINEARFULLGRID_HPP */
