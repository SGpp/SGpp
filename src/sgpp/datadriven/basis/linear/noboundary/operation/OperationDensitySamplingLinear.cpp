/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author A. Mo-Hellenbrand

#include "datadriven/basis/linear/noboundary/operation/OperationDensitySamplingLinear.hpp"
#include "datadriven/operation/OperationDensityConditional.hpp"
#include "datadriven/operation/OperationDensityMargTo1D.hpp"
#include "datadriven/operation/OperationDensitySampling1D.hpp"
#include "datadriven/operation/DatadrivenOpFactory.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg
{
namespace datadriven
{
	void OperationDensitySamplingLinear::doSampling(base::DataVector* alpha, base::DataMatrix* &samples, size_t num_samples) {

		size_t num_dims = this->grid->getStorage()->dim();

		//output matrix
		samples = new base::DataMatrix(num_samples, num_dims);

		size_t size = num_samples / num_dims;
		if (size <= 0) throw base::operation_exception("Error: # of dimensions greater than # of samples. Exit program...");

		for (size_t i=0; i < num_dims; i++) {

			// 1. marginalize to dim_i
			base::Grid* grid_i = NULL;
			base::DataVector* alpha_i = NULL;
			OperationDensityMargTo1D* marg1d = op_factory::createOperationDensityMargTo1D(*(this->grid));
			marg1d->margToDimX(alpha, grid_i, alpha_i, 0);
			delete marg1d;

			// 2. 1D sampling on dim_i
			base::DataVector* samples_i = new base::DataVector(1);
			OperationDensitySampling1D* samp1d = op_factory::createOperationDensitySampling1D(*grid_i);

			if (i == num_dims-1)
				size = num_samples - size*(num_dims-1);

			samp1d->doSampling1D(alpha_i, size, samples_i);
			delete samp1d;
			delete grid_i;
			delete alpha_i;

			// 3. for every sample on dim_i, do sampling
			base::DataVector* sampleVec = new base::DataVector(num_dims);
			for (size_t j=0; j < samples_i->getSize(); j++) {

				sampleVec->setAll(-1.0);
				sampleVec->set(i, samples_i->get(j));  //put samples_i[j] into i-th entry of samplesVec
				sampling_on_all_dims(this->grid, alpha, i, sampleVec);

				// copy result to output
				for (size_t k=0; k<num_dims; k++)
					samples->set(i*size + j, k, sampleVec->get(k));
			}

		}  // end for(i)

		return;
	}

	void OperationDensitySamplingLinear::sampling_on_all_dims(base::Grid* grid, base::DataVector* alpha, unsigned int dim_start, base::DataVector* &sampleVec) {

		unsigned int dim = sampleVec->getSize();
		unsigned int count = dim_start;

		if (dim_start == 0) {
			sampling_on_higher_dims(grid, alpha, dim_start, sampleVec, count);

		} else if (dim_start == dim - 1) {
			sampling_on_lower_dims(grid, alpha, dim_start, sampleVec, count);

		} else if ((dim_start > 0) && (dim_start < dim - 1)) {
			sampling_on_lower_dims(grid, alpha, dim_start, sampleVec, count);
			count = dim_start;
			sampling_on_higher_dims(grid, alpha, dim_start, sampleVec, count);

		} else {
			throw base::operation_exception("Error: dimension out of range. Exit program...");
		}

		return;
	}

	void OperationDensitySamplingLinear::sampling_on_lower_dims(base::Grid* g_in, base::DataVector* al_in, unsigned int dim_x, base::DataVector* &sampleVec, unsigned int &count) {

		// Recursive Call Stopper: when the lowest dimension (dim_0) is reached, stop
		if (count == 0) return;

		// Step 1: do conditional on dim_count
		base::Grid* g_out = NULL;
		base::DataVector* al_out = new base::DataVector(1);
		OperationDensityConditional* cond = op_factory::createOperationDensityConditional(*g_in);
		cond->doConditional(*al_in, g_out, *al_out, count, sampleVec->get(count));
		delete cond;

		count--;

		// Step 2: marginalize to dim_count for g_out
		base::Grid* g_1d = NULL;
		base::DataVector* al_1d = NULL;
		OperationDensityMargTo1D* marg1d = op_factory::createOperationDensityMargTo1D(*g_out);
		marg1d->margToDimX(al_out, g_1d, al_1d, count);
		delete marg1d;

		// Step 3: do sampling on dim_count
		base::DataVector* sample = new base::DataVector(1);
		OperationDensitySampling1D* samp = op_factory::createOperationDensitySampling1D(*g_1d);
		samp->doSampling1D(al_1d, 1, sample);
		delete samp;

		// Step 4: put the sample into output
		sampleVec->set(count, sample->get(0));
		delete sample;

		// repeat step 1-4 for the next dimension
		sampling_on_lower_dims(g_out, al_out, dim_x, sampleVec, count);
		delete g_out;
		delete al_out;

		return;
	}

	void OperationDensitySamplingLinear::sampling_on_higher_dims(base::Grid* g_in, base::DataVector* al_in, unsigned int dim_x, base::DataVector* &sampleVec, unsigned int &count) {

		// Recursive Call Stopper: when the lowest dimension (dim_0) is reached, stop
		if (count == (sampleVec->getSize() - 1)) return;

		// Step 1: do conditional on dim_count
		base::Grid* g_out = NULL;
		base::DataVector* al_out = new base::DataVector(1);
		OperationDensityConditional* cond = op_factory::createOperationDensityConditional(*g_in);
		cond->doConditional(*al_in, g_out, *al_out, dim_x, sampleVec->get(count));
		delete cond;

		count++;

		// Step 2: marginalize to dim_count for g_out
		base::Grid* g_1d = NULL;
		base::DataVector* al_1d = NULL;
		OperationDensityMargTo1D* marg1d = op_factory::createOperationDensityMargTo1D(*g_out);
		marg1d->margToDimX(al_out, g_1d, al_1d, dim_x);
		delete marg1d;

		// Step 3: do sampling on dim_count
		base::DataVector* sample = new base::DataVector(1);
		OperationDensitySampling1D* samp = op_factory::createOperationDensitySampling1D(*g_1d);
		samp->doSampling1D(al_1d, 1, sample);
		delete samp;

		// Step 4: put the sample into output
		sampleVec->set(count, sample->get(0));
		delete sample;

		// repeat step 1-4 for the next dimension
		sampling_on_higher_dims(g_out, al_out, dim_x, sampleVec, count);
		delete g_out;
		delete al_out;

		return;
	}
}
}








