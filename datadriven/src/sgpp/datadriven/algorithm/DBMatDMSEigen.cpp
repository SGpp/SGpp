// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatDMSEigen.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>

#include <iostream>

DBMatDMSEigen::DBMatDMSEigen() {}

DBMatDMSEigen::~DBMatDMSEigen() {}

void DBMatDMSEigen::solve(sgpp::base::DataMatrix& eigenVectors,
		          sgpp::base::DataVector& eigenValues,
			  sgpp::base::DataVector& alpha,
		          sgpp::base::DataVector& rhs, double lambda) {
  unsigned int n = eigenVectors.getNcols();
  //Create a matrix view for the eigenvectors
  gsl_matrix_view q = gsl_matrix_view_array(eigenVectors.getPointer(), n, n); 
  //Create a vector view for the right hand side
  gsl_vector_view b = gsl_vector_view_array(rhs.getPointer(), n); 
  //Create a vector view for the eigenvalues
  gsl_vector_view e = gsl_vector_view_array(eigenValues.getPointer(), n); 
  gsl_vector* res = gsl_vector_alloc(n);

  //Compute Q^T * b
  gsl_blas_dgemv(CblasTrans, 1., &q.matrix, &b.vector, 0., res);

  //Compute D^(-1) * Q^T * b (with D = E + lambda * I)
  for (unsigned int i = 0; i < n; i++) {
    gsl_vector_set(&e.vector, i,
		   1 / (gsl_vector_get(&e.vector, i) + lambda));
  }
  gsl_vector_mul(res, &e.vector);

  gsl_vector_view alphaView = gsl_vector_view_array(alpha.getPointer(), n);

  //Compute Q * D^(-1) * Q^T * b
  gsl_blas_dgemv(CblasNoTrans, 1., &q.matrix, res, 0., &alphaView.vector);
  gsl_vector_free(res);

}

#endif /* USE_GSL */


