#ifndef FULLGRIDSET_
#define FULLGRIDSET_
#include <iostream>
#include <vector>
#include "FullGrid.hpp"
#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "grid/GridStorage.hpp"
#include "grid/common/BoundingBox.hpp"

using namespace std;

namespace sg{
namespace combigrid {

/**
 * Represents a set of fullgrids, which can be combined into a sparse grid of given level and dimension
 * The fullgrids are the decomposition of a sparse grid/square root grid
 * Pointers to the fullgrids are stored in the vector grids, the fullgrids being ordered after first order norm in decreasing order
 * */
class FullGridSet{

public:

		typedef vector<FullGrid*> FULLGRID;

		typedef sg::base::GridStorage::index_type index_type;

		typedef index_type::index_type index_t;

		typedef index_type::level_type level_t;

	  /**Constructs a set of FullGrids for a sparsegrid of given dimension and level
	   * The initial values for the data vectors are 0
	   * @param dimension the dimension of the sparse grid
	   * @param n the level of the sparse grid*/
        FullGridSet(size_t dimension,size_t n, bool boundaryPoints = true)
        {
        	gridsWithBonudaryPoints_ = boundaryPoints;
			generate(dimension,n);
			dim=dimension;
			bbox=0;
        }

      /** Constructs a set of FullGrids for a sparsegrid of given dimension and level
       * The initial values for the data vectors are 0
       * @param dimension the dimension of the sparse grid
       * @param n the level of the sparse grid
       * @param type the type of the sparse grid(Linear,LinearBoundary, LinearTrapezoidBoundary,SquareRoot
       * @note don't use this constructor for super trapezoid sparse grids*/
        FullGridSet(size_t dimension,size_t n,const char *type , bool boundaryPoints = true)
        {
        	gridsWithBonudaryPoints_ = boundaryPoints;
        	generate(dimension,n,type);
        	dim=dimension;
        	bbox=0;
        }

       /** Constructs a set of FullGrids for a sparsegrid of given dimension and level for truncated trapezoid sparse grids
    	 * The initial values for the data vectors are 0
    	 * @param dimension the dimension of the sparse grid
    	 * @param n the level of the sparse grid
    	 * @param l_user the number of levels cut off*/
        FullGridSet(size_t dimension,size_t n,size_t l_user , bool boundaryPoints = true )
        {
        	gridsWithBonudaryPoints_ = boundaryPoints;
        	generate(dimension,n,l_user);
        	dim=dimension;
        	bbox=0;
        }

        FullGridSet(size_t dimension,size_t n,size_t *l_user , bool boundaryPoints = true )
     	{
        	gridsWithBonudaryPoints_ = boundaryPoints;
        	generate(dimension,n,l_user);
        	dim=dimension;
        	bbox=0;
        }

       /** Constructor for an adaptive set of fullgrids
         * @param dimension the dimension of the grids
         * @param levels the refinement level in every direction
         * @param type the type of the sparse grid the fullgrids add up to*/
        FullGridSet(size_t dimension, size_t *levels, const char* type , bool boundaryPoints = true )
        {
        	gridsWithBonudaryPoints_ = boundaryPoints;
        	generate(dimension,levels,type,0);
        	dim=dimension;
        	bbox=0;
        }

        /** Constructor for an adaptive set of fullgrids (for truncated trapezoid grids)
		 * @param dimension the dimension of the grids
		 * @param levels the refinement level in every direction
		 * @param type the type of the sparse grid the fullgrids add up to*/
        FullGridSet(size_t dimension, size_t *levels, size_t l_user , bool boundaryPoints = true )
        {
        	gridsWithBonudaryPoints_ = boundaryPoints;
        	std::vector<size_t> l_user_vect;
        	l_user_vect.resize(dimension,l_user);
        	generate(dimension,levels,&(l_user_vect[0]));
        	//generate(dimension,levels,"modlinearTrapezoidBoundary",l_user);
        	dim=dimension;
        	bbox=0;
        }

        /** Ctor for dimension adaptive T-CT, and also for dimension adaptive truncating
         * @param dimensions
         * @param levels level vector
         * @param l_user truncation vector*/
        FullGridSet(size_t dimension, size_t *levels, size_t* l_user , bool boundaryPoints = true )
		{
        	gridsWithBonudaryPoints_ = boundaryPoints;
        	generate(dimension,levels,l_user);
        	dim=dimension;
        	bbox=0;
		}

      /** Destructor of the fullgridset*/
      ~FullGridSet()
      {
           for (size_t i=0;i<size;i++)
               delete grids.at(i);
      }

      /** Lists all the coordinates of the FullGrids present in this set */
      void printGridSet()
      {
      	for (size_t i=0;i<size;i++)
      	{
      		cout<<i<<": ";
      		grids.at(i)->printCoords(dim);
      		cout<<"\n~~~~~~~~~~~~\n";
      	}
      }

	/** Initializes the data vectors of the FullGrids with the values from alpha, each value corresponding to the gridpoint with the same index from the sparsegrid
	 * Completes the fullgrids taking every gridpoint from the sparse grid
	 * @param storage pointer to the sparsegrid's storage which will be decomposed
	 * @param alpha the vector which contains the function values
	 * @note does the same thing as deCompose, but is slower */
	void initialize(sg::base::GridStorage *storage,sg::base::DataVector alpha);

	/** Prints the data values for each FullGrid */
	void printValues()
	{
		for (size_t i=0;i<size;i++){
			cout<<i<<": ";
		    grids.at(i)->printValues();
		}
	}

	/** Computes a result value baseds on the values stored in the value fields of the FullGrids
	  * Uses the combinational formula 4.20
	  * @return the linear combination of the values*/
      double combinedResult();

	 /** Returns a pointer to the set of FullGrids
	  * @return a pointer to the set of FullGrids*/
	  FULLGRID* getGrids()
      {
         return &grids;
      }

	 /** Returns the size of this set
	  * @return the total number of FullGrids*/
      size_t getSize()
      {
      	return size;
      }

       /** Allows acces to the n-th fullgrid
      *@param index the index of the fullgrid in the grids vector
      *@return the index-th fullgrid*/
      inline FullGrid& operator[](size_t index)
      {
	   return (*grids.at(index));
      }

     /**Allows acces to the n-th fullgrid
      *@param index the index of the fullgrid in the grids vector
      *@return a pointer to the index-th fullgrid*/
      FullGrid* at(unsigned int index)
      {
        return grids.at(index);
	  }

     /**Getter and setter for python's __getitem__ and __setitem__*/
     double get(size_t index,size_t pos)
	 {
	    return (*grids.at(index))[pos];
	 }

     /** set the value of a given full grid at a given position*/
     void set(size_t index,size_t pos, double val )
	 {
	    (*grids.at(index))[pos]=val;
	 }

     /**Returns the size of the n-th fullgrid
     *@param index the index of the fullgrid in the grids vector
     *@return the number of gridpoints in grids[index]*/
     size_t sizeOf(size_t index)
	 {
        return (*grids.at(index)).getSize();
	 }

     /** sets the combination coefficients*/
     void setCoefs(sg::base::DataVector &v)
     {
      	 for (size_t i=0;i<size;i++)
      	 {
       		 coefs[i]=v[i];
       	 }
     }

     /** returns the combination coefficients*/
     void getCoefs(sg::base::DataVector &v)
     {
       	 for (size_t i=0;i<size;i++)
       	 {
       		 v[i]=coefs[i];
       	 }
     }

	/** Decomposes a sparse grid into full grids, same functionality as the initialize method, but takes every point from the fullgrids and finds
	 * the corresponing sparse grid point.
	 * @param storage a pointer to the sparsegrid storage
	 * @param alpha the source datavector */
	void deCompose(sg::base::GridStorage *storage,sg::base::DataVector alpha);

	/** Combines a set of fullgrids into a sparse grid using the combination formula
	 * @param storage a pointer to the sparsegrid object
	 * @param alpha a pointer to the destination datavector */
	void reCompose(sg::base::GridStorage *storage,sg::base::DataVector *alpha);

	/** Evaluates all fullgrids and then combines the result
	 * @param p a datavector containing the coordinates of the point to evaluate the fullgrids in
	 * @return the value estimated on the set of fullgrids */
	double eval(sg::base::DataVector& p);

	/** Sets a bounding box for the fullgridset
	 * Is equivalent with a scaling of the axes
	 * @param boundingBox a pointer to the boudningbox of the set*/
	void setBoundingBox(sg::base::BoundingBox* boundingBox)
	{
		bbox=boundingBox;
		for (size_t i=0;i<size;i++)
			grids.at(i)->setBoundingBox(bbox);
	}

	/**  Gets the boundingbox of the fullgrids
	 * @return a pointer to the boundingbox*/
	sg::base::BoundingBox* getBoundingBox()
	{
		return bbox;
	}

	/** remove the full grids which have been added twice */
	void removeDuplicates();

protected:

	/** Generates all FullGrids for a given dimension and level(default type linearBoundary considered)
	 * They are computed by the formula (20)
	 * @param dim the dimension of the FullGrids
	 * @param n the level of the sparse grid for which the full grids are computed */
	void generate(size_t dim,size_t n);

    /** Generates all FullGrids for a given dimension and level and type
	 * They are computed by the formula (20)
	 * @param dim the dimension of the FullGrids
	 * @param n the level of the sparse grid for which the full grids are computed
     * @param type the type of the sparse grid which will be decomposed(returned by sg::base::Grid::getType())
     * @note currently implemented: linear, linearBoundary and linearTrapezoidBoundary*/
    void generate(size_t dim,size_t n,const char *type);

    /** Generate method for the SuperTrapezoid fullgrids
     * @param dim the dimension of the FullGrids
     * @param n the level of the sparse grid for which the full grids are computed
     * @param l_user the number of levels cut off*/
    void generate(size_t dim,size_t n,size_t l_user);

    /** the same function as above but with dimension adaptive truncating vector */
    void generate(size_t dim,size_t n,size_t* l_user);

     /** Generate method for adaptive grids
     * @param dim the dimension of the FullGrids
     * @param levels the refinement level in every direction
     * @param type the type of the sparse grid the grids add up to */
     void generate(size_t dim, size_t* levels, const char *type,size_t l_user);

     void generate(size_t dim, size_t* levels, size_t* l_user);

	/** Gets all FullGrids for which the individual levels sum up to a given number
	 * @param v auxiliar parameter, initially the empty vector, it collects the levels for every dimension
	 * @param dim the dimension of the full grids
	 * @param sum the remaining sum */
	void getsums(vector<level_t> *v,size_t dim,size_t sum);

	/** Same function for Linear Fullgrids without boundary */
    void getInnersums(vector<level_t> *v,size_t dim,size_t sum);

	/** Same function for trapezoid grids*/
    void getTrapezoidsums(vector<level_t> *v,size_t dim,size_t sum,size_t l_user);
    void getTrapezoidsums(vector<level_t> *v,size_t dim,size_t sum,size_t* l_user);

    /** Gets the fullgrid decomposition of a square root grid and stores it in the grids vector
     * @param v auxiliar parameter, initially empty vector
     * @param dim the dimension of the square root grid
     * @param maxlevel the level of the square root grid(h) */
     void getSquare(vector<level_t> *v, size_t dim,size_t maxlevel);

   /** The same methods for adaptive gridsets, with the ratio of every level relative to the coarsest level
    * We use the same method for all types of grids with boundary(except the squareRoot grid) */
   void getTrapezoidsums(vector<level_t> *v,size_t dim,int sum,int l_user,double* ratio);

   void getTrapezoidsums(vector<level_t> *v,size_t dim,int sum,size_t* l_user,double* ratio);

   void getInnersums(vector<level_t> *v,size_t dim,int sum,double* ratio);

   /** Computes the combination C(n,k)
   	 * @param n the first parameter for the combinational formula
   	 * @param k the second parameter for the combinational formula */
	static unsigned int combination(size_t n, size_t k)
	{
		if ((k==0)||(n==k)) return 1;
		else
		if ((k==1)||(n==k+1)) return n;
		else
		if ((k==2)||(n==k+2)) return n*(n-1)/2;
		else
		return combination(n-1,k)+combination(n-1,k-1);
	}


private:

       /**A vector of pointers to the fullgrids*/
	FULLGRID grids;

      /**The number of fullgrids in this set*/
	size_t size;

      /**The dimension of the fullgrids*/
	size_t dim;

     /**The type is 0 for linearBoundary fullgrids, 1 for linear fullgrids,2 for trapezoid fullgrids,3 for square root grids*/
    int gridType;

    /**The coefficient vector for the combi technique*/
    vector<double> coefs;

    /**The boundingbox of the grids*/
    sg::base::BoundingBox* bbox;

    /** flag to indicate if the created full grids should be with or wothout boundary points*/
    bool gridsWithBonudaryPoints_;
};
}
}
#endif /*FULLGRIDSET_*/
