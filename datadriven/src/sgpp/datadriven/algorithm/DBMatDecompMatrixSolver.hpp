// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DBMATDECOMPMATRIXSOLVER_HPP_
#define DBMATDECOMPMATRIXSOLVER_HPP_

#include <sgpp/solver/SGSolver.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

class DBMatDecompMatrixSolver: public sgpp::solver::SGSolver {
public:
	DBMatDecompMatrixSolver ();
	virtual ~DBMatDecompMatrixSolver();
};

#endif /* DBMATDECOMPMATRIXSOLVER_HPP_ */
