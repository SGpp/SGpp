/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#include "OperationMatrixLTwoDotPeriodic.hpp"
#include <sgpp/base/exception/data_exception.hpp>

namespace sg {
  namespace pde {

    OperationMatrixLTwoDotPeriodic::OperationMatrixLTwoDotPeriodic(sg::base::GridStorage* gridStorage) {
      this->gridStorage = gridStorage;
    }

    OperationMatrixLTwoDotPeriodic::~OperationMatrixLTwoDotPeriodic() {
      //TODO: FZ: delete gridStorage?
    }

    void OperationMatrixLTwoDotPeriodic::mult(sg::base::DataVector& alpha, sg::base::DataVector& result){
	  size_t nrows = gridStorage->size();
	  size_t ncols = gridStorage->size();

	  if (alpha.getSize() != ncols || result.getSize() != nrows) {
	    throw sg::base::data_exception("Dimensions do not match!");
	  }

	  size_t gridSize = gridStorage->size();
	  size_t gridDim = gridStorage->dim();

	  sg::base::DataMatrix level(gridSize, gridDim);
	  sg::base::DataMatrix index(gridSize, gridDim);

	  gridStorage->getLevelIndexArraysForEval(level, index);

	  sg::base::DataVector row(nrows);

	  for (size_t i = 0; i < gridSize; i++) {
		for (size_t j = 0; j < gridSize; j++) {
		  double res = 1;
		  for (size_t k = 0; k < gridDim; k++) {
			double lik = level.get(i, k);
			double ljk = level.get(j, k);
			double iik = index.get(i, k);
			double ijk = index.get(j, k);

			// i has always the lower level than j
			if(lik > ljk){
			  std::swap(lik,ljk);
			  std::swap(iik,ijk);
			}

			if(lik == 1){ // level 0
			  if(ljk > 2){
				lik = 2;
				iik = 1;

				if(ijk < ljk/2){
				  ijk = ljk/2 - ijk;
				}else{
				  ijk = ijk - ljk/2;
				}

				//Use formula for different overlapping ansatz functions:
				double diff = (ijk / ljk) - (iik / lik); // x_j_k - x_i_k
				double temp_res = fabs(diff - (1 / ljk))
								+ fabs(diff + (1 / ljk)) - fabs(diff);
				temp_res *= lik;
				temp_res = (1 - temp_res) / ljk;
				res *= temp_res;
			  }else {
				// if l2 == 0 => ljk = 1 =>
				// res = 1/3
				// if l2 == 1 => ljk = 2 =>
				// res = 1/6
				res *= 1.0 / (3 * ljk);
			  }
			}else if (lik == ljk) {
			  if (iik == ijk) { // case 4
				//Use formula for identical ansatz functions:
				res *= 2.0 / lik / 3;
			  } else { //case 0
				//Different index, but same level => ansatz functions do not overlap:
				res = 0.;
				break;
			  }
			} else {
			  if (std::max((iik - 1) / lik, (ijk - 1) / ljk)
				  >= std::min((iik + 1) / lik, (ijk + 1) / ljk)) {
				//Ansatz functions do not not overlap:
				res = 0.;
				break;
			  } else {
				//Use formula for different overlapping ansatz functions:

				double diff = (ijk / ljk) - (iik / lik); // x_j_k - x_i_k
				double temp_res = fabs(diff - (1 / ljk))
								+ fabs(diff + (1 / ljk)) - fabs(diff);
				temp_res *= lik;
				temp_res = (1 - temp_res) / ljk;
				res *= temp_res;
			  }
			}
		  }
		  row[j] = res;
		}
		//Standard matrix multiplication:
		double temp = 0.;
	    for (size_t j = 0; j < ncols; j++) {
		  temp += row[j] * alpha[j];
	    }
	    result[i] = temp;
	  }
	}
  }
}
