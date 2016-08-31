// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DBMATOFFLINE_H_
#define DBMATOFFLINE_H_

#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/base/grid/Grid.hpp>
//#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>

#ifdef USE_GSL
#include <gsl/gsl_permutation.h>
#endif /* USE_GSL */

using namespace std;
using namespace sgpp;

/**
 * Class that is used to decompose and store the left-hand-side matrix for the density based classification approach
 * (The classification is divided into two parts: the offline step that does not
 * depend on the actual data and the online step that depends on the data)
 */

class DBMatOffline {
protected:

	DBMatOffline();

	sgpp::datadriven::DBMatDensityConfiguration* config_; //configuration for this offline object
	sgpp::base::DataMatrix* lhsMatrix_; //stores the (decomposed) matrix
	bool constructed_; //If the matrix was built
	bool decomposed_; //If the matrix was decomposed
	bool ownsConfig_; //If the configuration has to be destroyed with the object
#ifdef USE_GSL
	gsl_permutation* perm_; //Stores the permutation that was applied on the matrix during decomposition
#endif /* USE_GSL */
	//An offline object either works on a full grid or on a hierarchical basis grid, so one of the following
	//pointers is non-null!
	sgpp::base::Grid* grid_;
	//combigrid::FullGrid<double>* fullgrid_;

	/**
	 * Method to initialize a sparse grid
	 */
	void InitializeGrid();
	/**
	 * Method to initialize a full grid
	 */
	//void InitializeFullGrid();

public:
	/**
	 * Constructor
	 *
	 * @param oc configuration for this offline object
	 */
	DBMatOffline(sgpp::datadriven::DBMatDensityConfiguration& oc);
	/**
	 * Constructor
	 *
	 * @param fname name of the file that stores the matrix + configuration
	 */
	DBMatOffline(string fname);

	/**
	 * Copy Constructor
	 * 
	 * The matrix needs to be already decomposed.
	 * 
	 * @param old Object to copy
	 */
	DBMatOffline(const DBMatOffline &old);

	/**
	 * Returns a pointer to the configuration
	 */
	sgpp::datadriven::DBMatDensityConfiguration* getConfig();

	/**
	 * Returns a pointer to the decomposed matrix
	 */
	sgpp::base::DataMatrix* getDecomposedMatrix();

	/**
	 * Returns a reference to the full grid (if this offline object uses a full grid, otherwise NULL)
	 */
	//combigrid::FullGrid<double>& getFullGrid();

	/**
	 * Returns a reference to the sparse grid (if this offline object uses a sparse grid, otherwise NULL)
	 */
	sgpp::base::Grid& getGrid();

	/**
	 * Returns a pointer to the sparse grid (if this offline object uses a sparse grid, otherwise NULL)
	 */
	sgpp::base::Grid* getGridPointer();

	/**
	 * Builds the right hand side matrix with or without the regularization term depending
	 * on the type of decomposition
	 */
	void buildMatrix();

	/**
	 * Decomposes the matrix according to the chosen decomposition type.
	 * The number of rows of the stored result depends on the decomposition type.
	 */
	void decomposeMatrix();

	/**
	 * Applies the permutation that might be done during the matrix decomposition to a vector
	 * (has to be applied to the right hand side vector if the system of equation should be solved with a decomposed matrix)
	 *
	 * @param b the vector that has to be permuted
	 */
	void permuteVector(sgpp::base::DataVector& b);
	
	/**
	 * Updates offline cholesky factorization based on coarsed (deletedPoints)
	 * and refined (newPoints) gridPoints
	 *
	 * @param deletedPoints list of indices of last coarsed points
	 * @param newPoints amount of refined points
	 */
	void choleskyModification(int newPoints,std::list<size_t> deletedPoints, double lambda);
	
	/**
	 * Updates the cholesky factor when a new grid point is added (e.g. refine)
	 *
	 * @param newCol DataVector with column to add to the system matrix
	 * @param size columns/rows of current Cholesky factor, necessary since the allocated memory is increased before the Cholesky factor is modified
	*/
	void choleskyAddPoint(sgpp::base::DataVector* newCol, size_t size);
	
	/**
	 * Permutes the rows of the cholesky factor based on permutations of the system matrix (e.g. coarsening)
	 *
	 * @param k "left" column to permutate
	 * @param l "right" column to permutate
	 * @param job = 2        => left circular shift 
	 *	  1,...,k-1,k,k+1, ..., l-1,l,l+1, ..,size  => 1,...,k-1,k+1, ..., l-1,l,k,l+1,..., size
	 * 	  job = 1       => right circular shift 
	 * 	  1,...,k-1,k,k+1, ..., l-1,l,l+1,...size  => 1,...,k-1,l,k,k+1, ..., l-1,l+1,...size
	*/
	void choleskyPermutation(int k, int l, int job);

	/**
	 * Prints the matrix onto standard output
	 */
	void printMatrix();

	/**
	 * Stores the matrix + configuration
	 *
	 * @param fname the file name
	 */
	void store(string fname);

	/**
	 * Loads matrix and configuration from a file
	 *
	 * @param fname the file name
	 */
	void load(string fname);

	/**
	 * Destructor
	 */
	virtual ~DBMatOffline();
private:
	/**
	 * Used for reading input files
	 */
	void Tokenize(string&, vector<string>&, string& delimiters);

};

#endif /* DBMATOFFLINE_H_ */
