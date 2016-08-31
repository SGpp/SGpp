// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DBMATDMSEigen_HPP_
#define DBMATDMSEigen_HPP_

#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>

/**
 * Class to solve the system of equations with a LU-decomposed matrix
 */
class DBMatDMSEigen: public DBMatDecompMatrixSolver {
public:
	/**
	 * (Empty) constructor
	 */
	DBMatDMSEigen();

	/**
	 * (Empty) destructor
	 */
	virtual ~DBMatDMSEigen();

	/**
	 * Solves a system of equations
	 *
	 * @param eigenVectors the eigendecomposed left hand side (the matrix contains the eigenvectors (rows 0...n) and eigenvalues (row n+1))
	 * @param alpha the vector of unknowns (the result is stored there)
	 * @param b the right hand vector of the equation system
	 */
	void solve(sgpp::base::DataMatrix& eigenVectors,
			sgpp::base::DataVector& eigenValues, sgpp::base::DataVector& alpha,
			sgpp::base::DataVector& rhs, double lambda);
};

#endif /* DBMATDMSEigen_HPP_ */
