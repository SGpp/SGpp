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
	size_t dims;
	SGPP::base::DataMatrix preparedDataset;
	STREAMING_OCL_INTERNAL_PRECISION *kernelDataset = nullptr;
	size_t datasetSize = 0;
	/// Member to store the sparse grid's levels for better vectorization
	STREAMING_OCL_INTERNAL_PRECISION* level = nullptr;
	/// Member to store the sparse grid's indices for better vectorization
	STREAMING_OCL_INTERNAL_PRECISION* index = nullptr;
	size_t gridSize = 0;
	/// Timer object to handle time measurements
	SGPP::base::SGppStopwatch myTimer;

	base::GridStorage* storage;

	float_t duration;

	OCLManager manager;
	OCLKernelImpl<STREAMING_OCL_INTERNAL_PRECISION> *kernel;
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
