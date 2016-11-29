// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <gsl/gsl_blas.h> 

#include <ctime>

DBMatOnlineDE::DBMatOnlineDE(double beta) :
  alpha_(nullptr), functionComputed_(false), b_save(nullptr),
  b_totalPoints(nullptr), beta_(beta), total_points(0),
  can_cv(false), test_mat(nullptr), test_mat_res(nullptr),
  lambda(0.), cv_logscale(false), normFactor(1.), o_dim(0) {
}

DBMatOnlineDE::~DBMatOnlineDE() {
  if (alpha_ != nullptr)
    delete alpha_;
  if (b_save != nullptr)
    delete b_save;
  if(b_totalPoints != nullptr)
    delete b_totalPoints;
}

void DBMatOnlineDE::readOffline(DBMatOffline* o) {
  offlineObject_ = o;
  functionComputed_ = false;
  if (b_save != nullptr)
    delete b_save;
  b_save = new sgpp::base::DataVector(offlineObject_->getDecomposedMatrix()->getNcols());
  b_totalPoints = 
    new sgpp::base::DataVector(offlineObject_->getDecomposedMatrix()->getNcols(), 0.0);
  lambda = offlineObject_->getConfig()->lambda_;
  o_dim = offlineObject_->getConfig()->grid_dim_;

  alpha_ = new sgpp::base::DataVector(offlineObject_->getDecomposedMatrix()->getNcols(), 0.0);
}

void DBMatOnlineDE::computeDensityFunction(sgpp::base::DataMatrix& m, bool save_b, bool do_cv,
					   std::list<size_t> *deletedPoints, 
                                           unsigned int newPoints) {
  if (m.getNrows() > 0) {
	
    sgpp::base::DataMatrix* lhsMatrix = offlineObject_->getDecomposedMatrix();

    //Compute right hand side of the equation:
    unsigned int numberOfPoints = m.getNrows();
    total_points++;
    sgpp::base::DataVector b(lhsMatrix->getNcols());
    b.setAll(0);
	
    std::unique_ptr<sgpp::base::OperationMultipleEval> B(
      sgpp::op_factory::createOperationMultipleEval(offlineObject_->getGrid(), m));
    sgpp::base::DataVector y(numberOfPoints);
    y.setAll(1.0);
    // Bt * 1
    B->multTranspose(y, b);
	
    //Perform permutation because of decomposition (LU)
    offlineObject_->permuteVector(b);

    //std::cout << b.getSize() << std::endl;

    if(save_b) {
      if(functionComputed_) {
	//double beta = std::max(beta_, (1./(double)total_points));
	//b.mult(beta);
	//b_save->mult(1.-beta);

	//Delete indices when grid got coarsend-> reduce 'b_save'
	if(deletedPoints != nullptr && !deletedPoints->empty()){
	  std::vector<size_t> v{std::begin(*deletedPoints), std::end(*deletedPoints)};
	  std::vector<size_t> v1(b_save->getSize() - deletedPoints->size());
	  size_t old_size = b_save->getSize();
		
	  size_t index_coarse = 0;
	  size_t index_remain = 0;
	  size_t temp;
	  for(size_t j = 0; j < old_size; j++){
	    temp = v[index_coarse];
	    if(temp == j){
	      index_coarse++;
	      continue;
	    }else{
	      v1[index_remain] = j;
	      index_remain++;
	    }
	  }
	  b_save->restructure(v1);
	  b_totalPoints->restructure(v1);	
	}
			
	//Expand 'b_save' when grid got refined
	if(newPoints > 0){
	  b_save->resizeZero(b.getSize());
	  b_totalPoints->resizeZero(b.getSize());
	}

	b.add(*b_save);
	//b.mult(beta);
      }
		
      //Update weighting based on processed data points
      for(size_t i = 0; i < b.getSize(); i++) {
	b_save->set(i, b.get(i));
	b_totalPoints->set(i, numberOfPoints + b_totalPoints->get(i));
	b.set(i, b_save->get(i) * (1. / b_totalPoints->get(i)));
      }
    }else {
      // 1 / M * Bt * 1
      b.mult(1. / numberOfPoints);
    }

    //Solve the system:
    if (alpha_ != nullptr)
      delete alpha_;
    alpha_ = new sgpp::base::DataVector(lhsMatrix->getNcols());

    //alpha_->resizeZero(lhsMatrix->getNcols());

    DBMatDecompostionType type = offlineObject_->getConfig()->decomp_type_;
    if (type == DBMatDecompLU) {
      DBMatDMSBackSub lusolver;
      lusolver.solve(*lhsMatrix, *alpha_, b);
    } else if (type == DBMatDecompEigen) {
      unsigned int n = lhsMatrix->getNcols();
      sgpp::base::DataVector e(n);
      lhsMatrix->getRow(n, e);
      DBMatDMSEigen esolver;
		
      if(can_cv && do_cv) {
	double best_crit = 0;
	double cur_lambda;
	for(int i = 0; i < lambda_step_; i++) {
	  cur_lambda = lambda_start_ + i*(lambda_end_ - lambda_start_)/(lambda_step_ - 1);
	  if(cv_logscale)
	    cur_lambda = exp(cur_lambda);
	  esolver.solve(*lhsMatrix, e, *alpha_, b, cur_lambda);
	  // double crit = computeL2Error();
	  double crit = resDensity(alpha_);
	  //std::cout << "cur_lambda: " << cur_lambda << ", crit: " << crit << std::endl;
	  if(i == 0 || crit < best_crit) {
	    best_crit = crit;
	    lambda = cur_lambda;
	  }

	}
      }
      esolver.solve(*lhsMatrix, e, *alpha_, b, lambda);
    } else if (type == DBMatDecompChol){
      DBMatDMSChol cholsolver;

      double old_lambda = lambda;
      //Perform cross-validation based on rank one up- and downdates
      //-> SHOULD NOT BE USED FOR LARGER GRID SETTINGS
      //ToDo: May be speed up by parallelization
      if(can_cv && do_cv) {
	double best_crit = 0;
	double cur_lambda;		
	for(int i = 0; i < lambda_step_; i++) {
	  cur_lambda = lambda_start_ + i*(lambda_end_ - lambda_start_)/(lambda_step_ - 1);
	  if(cv_logscale)
	    cur_lambda = exp(cur_lambda);
	  //std::cout << "Cur_lambda: " << cur_lambda << "  Old_lambda: " << old_lambda << std::endl;
	  //Solve for density declaring coefficients alpha based on changed lambda
	  cholsolver.solve(*lhsMatrix,*alpha_, b, old_lambda, cur_lambda);
	  old_lambda = cur_lambda;
	  double crit = resDensity(alpha_);
	  //std::cout << ", crit: " << crit << std::endl;
	  if(i == 0 || crit < best_crit) {
	    best_crit = crit;
	    lambda = cur_lambda;
	  }
	}
      }
      //Solve for density declaring coefficients alpha
      //std::cout << "lambda: " << lambda << std::endl;
      cholsolver.solve(*lhsMatrix, *alpha_, b, old_lambda, lambda);

    } else {
      throw sgpp::base::application_exception(
        "Unsupported decomposition type!");
    }
    functionComputed_ = true;
  }

}

double DBMatOnlineDE::resDensity(sgpp::base::DataVector* &alpha) {
  auto C = std::unique_ptr<sgpp::base::OperationMatrix>(
    sgpp::op_factory::createOperationIdentity(offlineObject_->getGrid()));
  sgpp::base::DataVector rhs(offlineObject_->getGrid().getSize());
  sgpp::base::DataVector res(offlineObject_->getGrid().getSize());
  sgpp::datadriven::DensitySystemMatrix SMatrix(offlineObject_->getGrid(), *test_mat, *C, 0.0);

  SMatrix.generateb(rhs);

  SMatrix.mult(*alpha, res);

  for(unsigned int i = 0; i < res.getSize(); i++)
    res[i] -= rhs[i];
  return res.l2Norm();
}

double DBMatOnlineDE::computeL2Error() {
  int nRows = test_mat_res->getNrows();
  sgpp::base::DataVector r(nRows);
  sgpp::base::DataVector tmp(test_mat->getNcols());
  for(int i = 0; i < nRows; i++) {
    test_mat->getRow(i, tmp);
    r[i] = this->eval(tmp, true);
  }
  double l2err = 0;
  for(int i = 0; i < nRows; i++) {
    l2err += (test_mat_res->get(i, 0) - r[i])*(test_mat_res->get(i, 0) - r[i]);
  }
  return sqrt(l2err)/nRows;
}

double DBMatOnlineDE::eval(sgpp::base::DataVector& p, bool force) {
  if (functionComputed_ || force == true) {
    double res;
    std::unique_ptr<sgpp::base::OperationEval> Eval(
      sgpp::op_factory::createOperationEval(offlineObject_->getGrid()));
    res = Eval->eval(*alpha_, p);
    return res*normFactor;
  } else {
    throw sgpp::base::data_exception("Density function not computed, yet!");
  }

}

sgpp::base::DataVector* DBMatOnlineDE::getAlpha() {
  return alpha_;
}

void DBMatOnlineDE::updateAlpha(std::list<size_t>* deletedPoints, unsigned int newPoints) {
  if (alpha_ != nullptr && deletedPoints != nullptr && !deletedPoints->empty()) {  
    std::vector<size_t> deletedPoints_{std::begin(*deletedPoints), std::end(*deletedPoints)};
    sgpp::base::DataVector* newAlpha =
      new sgpp::base::DataVector(alpha_->getSize()-deletedPoints->size()+newPoints);
    for (size_t i = 0; i < alpha_->getSize(); i++) {
      if (std::find(deletedPoints_.begin(), deletedPoints_.end(), i) != deletedPoints_.end()) {
        continue;
      }
      else {
        newAlpha->append(alpha_->get(i));  
      }      
    }    
    // delete old alpha
    delete alpha_;
    // set new alpha
    alpha_ = newAlpha;
  }
  if (newPoints > 0) {
    alpha_->resizeZero(alpha_->getSize()+newPoints);
  }
}

bool DBMatOnlineDE::isComputed() {
  return functionComputed_;
}

void DBMatOnlineDE::setCrossValidationParameters(int lambda_step, double lambda_start,
						 double lambda_end, sgpp::base::DataMatrix *test,
						 sgpp::base::DataMatrix *test_cc, bool logscale) {
  lambda_step_ = lambda_step;
  cv_logscale = logscale;
  if(cv_logscale) {
    lambda_start_ = std::log(lambda_start);
    lambda_end_ = std::log(lambda_end);
  }
  else {
    lambda_start_ = lambda_start;
    lambda_end_ = lambda_end;
  }
  if(test != nullptr)
    test_mat = test;
  if(test_cc != nullptr)
    test_mat_res = test_cc;
  can_cv = true;
}

double DBMatOnlineDE::getBestLambda() {
  return lambda;
}

void DBMatOnlineDE::setBeta(double beta) {
  beta_ = beta;
  return;
}

double DBMatOnlineDE::getBeta() {
  return beta_;
}

double DBMatOnlineDE::normalize(size_t samples) {
  this->normFactor = 1.;
  double sum = 0.;
  sgpp::base::DataVector p(this->o_dim);
  srand(time(nullptr));
  for(size_t i = 0; i < samples; i++) {
    for(size_t j = 0; j < this->o_dim; j++)
      p[j] = (static_cast<double>(rand())/RAND_MAX);
    sum += this->eval(p);
  }
  return this->normFactor = samples/sum;
}

#endif /* USE_GSL */
