#ifndef FULLGRID_
#define FULLGRID_
#include <iostream>
#include <vector>
#include "math.h"
#include "grid/GridStorage.hpp"
#include "grid/common/BoundingBox.hpp"
using namespace std;
namespace sg{

/** the arrays which contains the precalculated values of power of two*/
const int powOfTwo[30]= {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216,33554432,67108864,134217728,268435456};
const double negPowOfTwo[30]= {1.0/1,1.0/2,1.0/4,1.0/8,1.0/16,1.0/32,1.0/64,1.0/128,1.0/256,1.0/512,1.0/1024,1.0/2048,1.0/4096,1.0/8192,1.0/16384,1.0/32768,1.0/65536,1.0/131072,1.0/262144,1.0/524288,1.0/1048576,1.0/2097152,1.0/4194304,1.0/8388608,1.0/16777216,1.0/33554432,1.0/67108864,1.0/134217728,1.0/268435456};

/**
 * Represents a single full grid of given dimension and level
 * The index of a gridpoint in the full grid is given by the formula ind=i0+i1*N0+i2*N0*N1+...+id*N0*N1*N2*...*Nd, where i0,i1,i2,... are the indexes in every dimension, and Nk is the number of gridpoints in direction k
 * Nk=2^level[k]+1 for every k for the LinearBoundaryFullGrid(default) and Nk=2^level[k]-1 for the LinearGrid without boundary
 * */
class FullGrid{

	
  public:

	  // to simplify the types we declare here the types
	  typedef GridStorage::index_type index_type;
	  typedef index_type::index_type index_t;
 	  typedef index_type::level_type level_t;

    /** Factory method for creation of a linear fullgrid without boundary
    *@param dim the dimension of the grid
    *@param inlevel a pointer to a vector which contains the levels of the fullgrid
    *@return a pointer to a linear fullgrid*/
    static FullGrid* createLinearFullGrid(size_t dim, vector<level_t> *inlevel);    

    /** Factory method for creation of a linear fullgrid with boundary
    *@param dim the dimension of the fullgrid
    *@param inlevel a pointer to a vector which contains the levels of the fullgrid
    *@return a pointer to a linear boundary fullgrid */
    static FullGrid* createLinearBoundaryFullGrid(size_t dim, vector<level_t> *inlevel);

    /** Copy constructor for the fullgrid
     * @param g the fullgrid to copy */
    FullGrid(const FullGrid &g):dim(g.dim),size(g.size),vec(g.vec),boundingBox(g.boundingBox)
	{	
		level=new level_t[dim];
		for (size_t i=0;i<dim;i++)
		{
			level[i]=g.level[i];
		}			
    	 intersect=new double[2*dim];
    	 aindex=new int[dim];
	}

	/** Creates a fullgrid skeleton of given dimension and level
	 * @param indim the dimension of the grid*/
    FullGrid(size_t indim): dim(indim),size(1),boundingBox(0)
    {          
      level=new level_t[dim];
 	  intersect=new double[2*dim];
 	  aindex=new int[dim];
    }

    /** Ctor for dimension adaptive full grid
     * @param */
    FullGrid(size_t indim, vector<level_t> *inlevel): dim(indim),size(1),boundingBox(0){
       level=new level_t[dim];
       for (size_t i=0;i<dim;i++){
			  level[i]=inlevel->at(i);
			  size=size*(powOfTwo[inlevel->at(i)]+1);
		}         		
       vec.resize(size,0);
       intersect=new double[2*dim];
       aindex=new int[dim];
     }

    /** Destructor*/
    ~FullGrid()
    {    	
      delete[] level;
      delete[] intersect;
      delete[] aindex;
    }

    /** Copy constructor */
    FullGrid operator=(const FullGrid &fg)
    {
      if (this==&fg) return *this;
      if (dim!=fg.dim)
          {
	      dim=fg.dim;
	      level=new level_t[dim];
	      size=fg.size;
          vec.resize(size,0);
	  }       
      for (size_t i=0;i<dim;i++)		
			level[i]=fg.level[i];		
      for (size_t i=0;i<size;i++)
		        vec.at(i)=fg.vec.at(i);
 	 intersect=new double[2*dim];
 	 aindex=new int[dim];
 	 boundingBox=fg.boundingBox;
 	 return *this;
    }

    /** Returns the number of gridpoints
     * @return the size of the grid*/
    size_t getSize(){ return size;}

    /** Returns a pointer to the level vector
     * @return the level of the full grid */
    inline level_t *getLevel(){return level;}

    inline level_t getLevel(size_t index)
    {
    	return level[index];
    }

    /** Returns the dimension of the grid
     * @return the dimension */
    size_t getDim(){return dim;}    

    /** Allows acces to the i-th gridpoint value by index in the vector
     * @param index the index of the gridpoint
     * @return acces to the gridpoint value at the position index */
	inline double& operator[](size_t index)
    {
       return vec.at(index);
    }
	/** Allows acces to the gridpoints by it's indexes in every dimension
	 * @param index the index of the gridpoint in every dimension
	 * @return acces to the gridpoint value at the given coordinate */
    virtual inline double& at(size_t *index){
    	size_t ind=0;
		for (int i=dim-1;i>=0;i--)
    	    	ind=(powOfTwo[level[i]]+1)*ind+index[i];
		return vec.at(ind);
    }

     /** Gives the coordinates of the gridpoint in a cartesian coord system and puts them in the corresponding Datavector
     * @param index the index of the gridpoint in the vector
     * @param v a DataVector-will contain double values from the interval [0,1] representing the coordinates in all directions */
    virtual void getCoords(size_t index,DataVector &v);

    /** Gives the coordinates of the gridpoint in a cartesian coord system
     * @param d the dimension for which we search the coordinate
     * @param index the index of the gridpoint in the vector
     * @return a double value from the interval [0,1] representing the coordinate in direction d */
    virtual inline double getCoord(size_t d,size_t index){
    	int ind=index;
    	double c=0;

    	/** Computes coordinate of gridpoint in direction d as: [[index/(N0*N1*N2...*Nd-1)] %Nd] *2^(-level[d])
    	 *  Ni=2^level[i]+1*/
    	for (size_t i=0;i<d;i++)
    	{
    		 ind=ind / (powOfTwo[level[i]]+1);
    	}
        ind=ind%((int)powOfTwo[level[d]]+1);
        if (ind==0)
			{ c=0.0; }
        else
        	{ c=((double)ind)*negPowOfTwo[level[d]]; }
        if (boundingBox!=0)
	    {
		    c=c*boundingBox->getIntervalWidth(d)+boundingBox->getIntervalOffset(d);
	    }
        return c;
    }

    /** Generates a string with all coordinates of the grid point.
     * The accuracy is up to 6 digits, i.e. beginning with level 8 there are rounding errors.
     * @return returns a string with the coordinates of the grid point separated by whitespace*/
    virtual std::string getCoordsString(size_t index);

    /** Prints all gridpoint values */
    void printValues()
    {
    	for (size_t i=0;i<size;i++)
    	{    		
    		cout<<vec.at(i);
    		cout<<" ";
    	}
    	cout<<"\n";
    }

    /** Prints the level of this fullgrid */
    void printCoords(size_t dim)
    {
       cout<<getType()<<": ";
       for (size_t i=0;i<dim;i++)
    	{
    		cout<<level[i]<<" ";    		
    	}    	
        cout<<"\n";    	
    }

    /**Returns the 1-norm of the level vector
     * @return the sum of levels */
    virtual size_t norm1()
    {
    	size_t sum=0;
    	for (size_t i=0;i<dim;i++)
    		sum+=level[i];
    	return sum;
    }

    /** Returns the value assigned to this fullgrid after evaluation */
    double &val()
    {
    	return value;
    }

    /**Initializes the double data vector */
    double get(size_t index)
    {
          return vec.at(index);
    }

    /** Set and get for python's __getitem__ and __setitem__*/
    void set(size_t index,double value)
    {
          vec.at(index)=value;
    }

    /**Initializes the values of the data vector from a given vector
     *@param v the input vector*/
     void fill(DataVector &v)
     {
     	for (size_t i=0;i<vec.size();i++)
     		vec.at(i)=v[i];
     }

    /** Returns the type of the grid-default linearBoundary */
    virtual const char* getType()
    {
      return "linearBoundary";
    }

    /** Returns the number of gridpoints in direction d  */
    virtual inline size_t length(size_t d)
    {
       return powOfTwo[level[d]]+1;
    }

    /** Returns the index of the first gridpoint on direction d */
    virtual inline size_t startindex()
    {
       return 0;
    }

     /** Evaluates the value of the function in a given point
      * @param p the coordinates of the point
      * @return the approximation of the function in the point p */
     virtual double eval(DataVector& p);

 	/** Sets a bounding box for the fullgrid
 	 * Is equivalent with a scaling of the axes
 	 * @param boundingBox a pointer to the boudningbox of the grid
 	 * */
     void setBoundingBox(BoundingBox* bBox)
     {
    	 boundingBox=bBox;
     }

	/** Gets the boundingbox of the fullgrids
	 * @return a pointer to the boundingbox*/
     BoundingBox* getBoundingBox()
	 {
    	 return boundingBox;
	 }

     /** Two fullgrids are equal if their levels equal on every dimension
      * Obs. It is not necessary to have the same values at the gridpoints
      * @param fg the fullgrid to compare with
      * @return true if all the levels of this grid are equal to the levels of fg*/
     inline bool equals(FullGrid& fg)
     {
	  for (size_t i=0;i<dim;i++)
		  if (fg.level[i]!=level[i]) return false;
	  return true;
     }


  protected:

  /**the dimension of the FullGrid*/
    size_t dim;

  /**the number of gridpoints*/
    size_t size;

  /**the double data vector containing the function values for every gridpoint*/
    vector<double> vec;

  /**the levels of the fullgrid in an array of dimension dim*/
    level_t *level;

  /**the function value for the FullGrid in the given point*/
    double value;

  /**The boundingbox of the grid*/
   BoundingBox* boundingBox;

  /**Auxiliar variables for use with eval*/
	 double *intersect;
	 int *aindex;

};}

#endif /*FULLGRID_*/
