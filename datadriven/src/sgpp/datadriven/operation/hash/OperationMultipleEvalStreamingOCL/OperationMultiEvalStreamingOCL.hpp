// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include "OCLKernelImpl.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class OperationMultiEvalStreamingOCL: public base::OperationMultipleEval {
protected:
	SGPP::base::DataMatrix preparedDataset;
	/// Member to store the sparse grid's levels for better vectorization
	SGPP::base::DataMatrix* level_ = nullptr;
	/// Member to store the sparse grid's indices for better vectorization
	SGPP::base::DataMatrix* index_ = nullptr;
	/// Timer object to handle time measurements
	SGPP::base::SGppStopwatch myTimer_;

	base::GridStorage* storage;

	float_t duration;

	OCLManager manager;
	OCLKernelImpl *kernel;
public:

	OperationMultiEvalStreamingOCL(base::Grid& grid, base::DataMatrix& dataset);

	~OperationMultiEvalStreamingOCL();

	void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result)
			override;

	void multTranspose(SGPP::base::DataVector& source,
	SGPP::base::DataVector& result) override;

	float_t getDuration();

	void prepare() override;

private:
	size_t padDataset(SGPP::base::DataMatrix& dataset);

	void recalculateLevelAndIndex();
};

}
}
