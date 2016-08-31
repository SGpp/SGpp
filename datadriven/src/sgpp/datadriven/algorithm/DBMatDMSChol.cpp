// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>
#endif /* USE_GSL */

#include <math.h>
#include <iostream>
#include <ctime>

DBMatDMSChol::DBMatDMSChol() {
}

DBMatDMSChol::~DBMatDMSChol() {
}

void DBMatDMSChol::solve(sgpp::base::DataMatrix& DecompMatrix,
		sgpp::base::DataVector& alpha, sgpp::base::DataVector& b, double lambda_old, double lambda_new) {


	unsigned int size  = DecompMatrix.getNcols(); 
	//Performe Update based on Cholesky - afterwards perform n (GridPoints) many rank-One-updates

	double lambda_up = lambda_new - lambda_old;

	//If regularization paramter is changed enter
	if(lambda_up != 0){
		sgpp::base::DataVector* lambdaModification = new sgpp::base::DataVector(size, 0.0);
		//In case lambda is increased apply Cholesky rank one updates
		if(lambda_up > 0){
			for(int i = 0; i < size; i++){
				lambdaModification->set(i, sqrt(lambda_up));
				choleskyUpdate(DecompMatrix, lambdaModification, true);
				lambdaModification->setAll(0.0);
			}
		}else if(lambda_up < 0){
		//In case lambda is decreased apply Cholesky rank one downdates
			for(int i = 0; i < size; i++){
				lambdaModification->set(i, sqrt(fabs(lambda_up)));
				choleskyDowndate(DecompMatrix, lambdaModification, true);
				lambdaModification->setAll(0.0);
			} 
		}

		delete lambdaModification;
	}
	
	//Solve (R + lambda * I)alpha = b to obtain density declaring coefficents alpha.	

	//Forward Substitution:
	sgpp::base::DataVector y(size);
	for (unsigned int i = 0; i < size; i++) {
		y[i] = b[i];
		for (unsigned int j = 0; j < i; j++) {
			y[i] -= DecompMatrix.get(i, j) * y[j];
		}
		y[i] /= DecompMatrix.get(i, i);
	}

	//Backward Substitution:
	for (int i = size - 1; i >= 0; i--) {
		alpha[i] = y[i];
		for (unsigned int j = i + 1; j < size; j++) {
			alpha[i] -= DecompMatrix.get(j, i) * alpha[j];
		}
		alpha[i] /= DecompMatrix.get(i, i);
	}

	return;
}

// Implement cholesky Update for given Decomposition and update vector
void DBMatDMSChol::choleskyUpdate(sgpp::base::DataMatrix& DecompMatrix, sgpp::base::DataVector* update, bool do_cv){
#ifdef USE_GSL
	int i;

	unsigned int size = DecompMatrix.getNrows();

	if(update->getSize() != size){
		throw sgpp::base::data_exception("choleskyUpdate::Size of DecomposedMatrix and updateVector don´t match...");
	}
	
	gsl_matrix_view m = gsl_matrix_view_array(DecompMatrix.getPointer(), size, size); //Create GSL matrix view for update procedures
	gsl_vector_const_view vvec = gsl_vector_const_view_array(update->getPointer(), update->getSize());

	//Define and declare Workingvector, Cosine- and Sinevector
	gsl_vector*  wkvec = gsl_vector_calloc(update->getSize());
	gsl_vector*  svec  = gsl_vector_calloc(update->getSize());
	gsl_vector*  cvec  = gsl_vector_calloc(update->getSize());


	//Generate Givens rotations, update L 
  	//Copy Values of updateVector into WorkingVector
	gsl_blas_dcopy(&vvec.vector, wkvec);
	
	double temp;
	double* tbuff = m.matrix.data;
	bool first_notZero = true;
	for (size_t i = 0; i < size-1; i++){
		if (*tbuff == 0.0 && wkvec->data[i] == 0.0) {
			throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite");
    		}else if (wkvec->data[i] == 0.0 && do_cv == true){
		//If cross validation is applied the first (n-1) entries in the n-th step are zero and therefore can be skipped.
			tbuff += (size+1);
			continue;
		}else if (wkvec->data[i] == 0.0 && first_notZero == true){
			first_notZero = false;
			tbuff += (size+1);
			continue;
		}
		//Determine givens rotation
		gsl_blas_drotg(tbuff, wkvec->data+i, cvec->data+i, svec->data+i);
		if ((temp = *tbuff) < 0.0) {
      			*tbuff=-temp;
			cvec->data[i]=-cvec->data[i];
			svec->data[i]=-svec->data[i];
    		} else if (temp == 0.0) {
			throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite");
    		}
		
    		//Access columns to modify via Givens rotations
		gsl_vector_view mat_sub  = gsl_matrix_subcolumn(&m.matrix, i, i + 1, size - i - 1);
		gsl_vector_view wkvec_sub = gsl_vector_subvector(wkvec, i + 1, size - i - 1);

		//Allpy Givens rotation to mat_sub und wkvec_sub
		for(size_t j = 0; j < size - i - 1; j++){
			double x = mat_sub.vector.data[mat_sub.vector.stride * j];
			double y = wkvec_sub.vector.data[j];
			mat_sub.vector.data[mat_sub.vector.stride * j] = cvec->data[i] * x + svec->data[i] * y;
			wkvec_sub.vector.data[j] = -svec->data[i] * x + cvec->data[i] * y;
		}
   		tbuff += (size+1);
	}

	//Apply changes to N-th (last) diagonal element
	//Is outsourced, since only the diagonal element is modified.
	if (*tbuff != 0.0 || wkvec->data[size-1] != 0.0) {
		gsl_blas_drotg(tbuff, wkvec->data+(size-1), cvec->data+i, svec->data+i);
    		if ((temp = *tbuff) < 0.0) {
      			*tbuff=-temp; cvec->data[i]=-cvec->data[i]; svec->data[i]=-svec->data[i];
    		} else if (temp == 0.0){
			 throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite"); 
		}
  	} else {
		throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite"); 
	}
	
	gsl_vector_free(wkvec);
	gsl_vector_free(svec);
	gsl_vector_free(cvec);

	return; 
#endif /* USE_GSL */
}

// Implement cholesky Downdate for given Decomposition and update vector
void DBMatDMSChol::choleskyDowndate(sgpp::base::DataMatrix& DecompMatrix, sgpp::base::DataVector* downdate, bool do_cv){
#ifdef USE_GSL
	unsigned int size = DecompMatrix.getNrows();

	if(downdate->getSize() != size){
		throw sgpp::base::data_exception("choleskyDowndate::Size of DecomposedMatrix and updateVector don´t match...");
	}
	
	gsl_matrix_view m = gsl_matrix_view_array(DecompMatrix.getPointer(), size, size); //Create GSL matrix view for update procedures
	gsl_vector_view vvec = gsl_vector_view_array(downdate->getPointer(), downdate->getSize());

	//Define and declare Workingvector, Cosine- and Sinevector
	gsl_vector*  wkvec = gsl_vector_calloc(downdate->getSize());
	gsl_vector*  svec  = gsl_vector_calloc(downdate->getSize());
	gsl_vector*  cvec  = gsl_vector_calloc(downdate->getSize());
	
	//Compute p (if not given) 
  	//Copy Values of updateVector into WorkingVector
	gsl_blas_dcopy(&vvec.vector, wkvec);
	
    	//Solve La = downdate and save a in vvec.vector
	gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &m.matrix, &vvec.vector);
	
	//Generate Givens rotations
	double rho;
	gsl_blas_ddot (&vvec.vector, &vvec.vector, &rho);
	rho = 1 - rho;
	
	//Represents first index of downdate vector with entrie unequal to zero
	int cv_first_zero = 0;
	 if (rho<=0.0){
    		throw sgpp::base::data_exception("choleskyDowndate::Matrix not numerical positive definite"); 
  	} else {
    		rho=sqrt(rho);
    			for (int i = size-1; i>=0; i--) {
				//If this method is applied in cross-validation do_cv == true and shortcuts can be used
				if(vvec.vector.data[i] == 0 && do_cv == true){
					cv_first_zero = i + 1;
					break;
				}
		      		//Determine Givens rotation 
				gsl_blas_drotg(&rho, vvec.vector.data+i, cvec->data+i, svec->data+i);
		      		//rho must remain positive
		      		if (rho < 0.0) {
					rho=-rho; cvec->data[i]=-cvec->data[i]; svec->data[i]=-svec->data[i];
		      		}
			}
    	}

	//rho should be 1 now
	gsl_vector_set_zero(wkvec);
	
	double* tbuff = m.matrix.data +((size-1)*(size+1));

	//Apply calculated Givens rotations to current Cholesky factor
	for (int i=size-1; i>=cv_first_zero; i--) {
		if (*tbuff<=0.0) {
			throw sgpp::base::data_exception("choleskyDowndate::Matrix not numerical positive definite"); 
                }
		//Access corresponding columns
		gsl_vector_view mat_sub  = gsl_matrix_subcolumn(&m.matrix, i, i, size - i);
		gsl_vector_view wkvec_sub = gsl_vector_subvector(wkvec, i, size - i);

		//Apply Givens rotation to mat_sub and wkvec_sub
		gsl_blas_drot(&wkvec_sub.vector, &mat_sub.vector, cvec->data[i], svec->data[i]);
		double diag = gsl_matrix_get (&m.matrix, i, i);
		//Ensure diagonal stays positive
		if (diag<0.0) {
			rho=-1.0;
			gsl_matrix_set(&m.matrix, i, i, rho * diag);
      	        } else if (diag == 0.0) {
			throw sgpp::base::data_exception("choleskyDowndate::Matrix not numerical positive definite"); 
      	        }
	      tbuff-=(size+1); 
	}
	
	//Workingvector should equals v now

	gsl_vector_free(wkvec);
	gsl_vector_free(svec);
	gsl_vector_free(cvec);

	return;
#endif /* USE_GSL */	
}

