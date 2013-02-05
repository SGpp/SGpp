/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author A. Mo-Hellenbrand

#include "datadriven/basis/linear/noboundary/operation/OperationDensityMargTo1DLinear.hpp"
#include "datadriven/operation/OperationDensityMarginalize.hpp"
#include "datadriven/operation/DatadrivenOpFactory.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg
{
namespace datadriven
{
  void OperationDensityMargTo1DLinear::margToDimX(base::DataVector* alpha, base::Grid* &grid_x, base::DataVector* &alpha_x, unsigned int dim_x) {

	unsigned int dim_h = this->grid->getStorage()->dim() - 1;

	if (dim_x == 0) {
		margToDim0(this->grid, alpha, grid_x, alpha_x);

	} else if (dim_x == dim_h) {
		margRemoveLowerDims(this->grid, alpha, grid_x, alpha_x, dim_x, dim_h);

	} else if ((dim_x > 0) && (dim_x < dim_h)) {
		base::Grid* grid_tmp = NULL;
		base::DataVector* alpha_tmp = new base::DataVector(1);

		margRemoveLowerDims(this->grid, alpha, grid_tmp, alpha_tmp, dim_x, dim_h);
		margToDim0(grid_tmp, alpha_tmp, grid_x, alpha_x);

		delete grid_tmp;
		delete alpha_tmp;

	} else {
		throw base::operation_exception("Error: dimension out of range. Exit program...");
	}

	return;
  }

  void OperationDensityMargTo1DLinear::margRemoveLowerDims(base::Grid* g_in, base::DataVector* al_in, base::Grid* &g_out, base::DataVector* &al_out, unsigned int dim_x, unsigned int dim_h) {

	unsigned int target_dims = dim_h + 1 - dim_x;

	if (g_in->getStorage()->dim() == target_dims) {
		g_out = g_in;
		al_out = al_in;
		return;
	}

	// marginalize: remove one lower dim
	OperationDensityMarginalize* marg = op_factory::createOperationDensityMarginalize(*g_in);
	base::Grid* g_tmp = NULL;
	base::DataVector* al_tmp = new base::DataVector(1);

	marg->doMarginalize(*al_in, g_tmp, *al_tmp, 0);
	delete marg;

	// check the new grid's dimensions
	if (g_tmp->getStorage()->dim() == target_dims) {
		g_out = g_tmp;
		al_out = al_tmp;
	} else {
		margRemoveLowerDims(g_tmp, al_tmp, g_out, al_out, dim_x, dim_h);
		delete g_tmp;
		delete al_tmp;
	}

	return;
  }

  void OperationDensityMargTo1DLinear::margToDim0(base::Grid* g_in, base::DataVector* al_in, base::Grid* &g_out, base::DataVector* &al_out) {

	if (g_in->getStorage()->dim() == 1) {
		g_out = g_in;
		al_out = al_in;
		return;
	}

	// marginalize: remove dim 1
	OperationDensityMarginalize* marg = op_factory::createOperationDensityMarginalize(*g_in);
	base::Grid* g_tmp = NULL;
	base::DataVector* al_tmp = new base::DataVector(1);

	marg->doMarginalize(*al_in, g_tmp, *al_tmp, 1);
	delete marg;

	// check the new grid's dimensions
	if (g_tmp->getStorage()->dim() == 1) { //if dim = 1, return
		g_out = g_tmp;
		al_out = al_tmp;
	} else {                               //if dim > 1, marginalize further
		margToDim0(g_tmp, al_tmp, g_out, al_out);
		delete g_tmp;
		delete al_tmp;
	}

	return;
  }

}
}




