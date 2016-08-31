// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DBMatDMSBackSub_HPP_
#define DBMatDMSBackSub_HPP_

#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>

/**
 * Class to solve the system of equations with a LU-decomposed matrix
 */

class DBMatDMSBackSub: public DBMatDecompMatrixSolver {
public:
	/**
	 * (Empty) constructor
	 */
	DBMatDMSBackSub();

	/**
	 * (Empty) destructor
	 */
	virtual ~DBMatDMSBackSub();

	/**
	 * Solves a system of equations
	 *
	 * @param DecompMatrix the LU decomposed left hand side
	 * @param alpha the vector of unknowns (the result is stored there)
	 * @param b the right hand vector of the equation system
	 */
	void solve(sgpp::base::DataMatrix& DecompMatrix, sgpp::base::DataVector& alpha,
				sgpp::base::DataVector& b);
};
#endif /* DBMatDMSBackSub_HPP_ */
