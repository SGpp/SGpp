// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/PrimalDualSVM.hpp>
#include <sgpp/datadriven/operation/hash/OperationEvalSGKernel/OperationEvalSGKernel.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace datadriven {


PrimalDualSVM::PrimalDualSVM(size_t dim, size_t dataDim, size_t budget, bool useBias) 
    : svs(new sgpp::base::DataMatrix(0,dataDim)),       
      alphas(new sgpp::base::DataVector(0)),
      norms(new sgpp::base::DataVector(0)), 
      w(new sgpp::base::DataVector(dim, 0.0)),    
      w2(new sgpp::base::DataVector(dim, 0.0)), 
      bias(0.0),
      budget(budget),
      useBias(useBias) {

  svs->setInc(budget);
}

PrimalDualSVM::~PrimalDualSVM() {}


double PrimalDualSVM::predictRaw(sgpp::base::Grid& grid, sgpp::base::DataVector& x, 
                                 size_t dataDim, bool trans) {
  //SGKernel evaluation operation
  std::unique_ptr<datadriven::OperationEvalSGKernel> opSGKernel =
    op_factory::createOperationEvalSGKernel(grid);       
  
  sgpp::base::DataVector xTrans(grid.getSize()); 
  if (trans) {
    xTrans = x;
  }
  else {
    opSGKernel->phi(x, xTrans, dataDim);
  }
  double res = w->dotProduct(xTrans);
  if (useBias) {
    res += bias;
  }
  return res;
}


int PrimalDualSVM::predict(sgpp::base::Grid& grid, sgpp::base::DataVector& x, 
                           size_t dataDim) {
  bool sign = std::signbit(this->predictRaw(grid, x, dataDim)); 
  if (sign) {
    return -1; 	  
  }
  else {
    return 1;  
  }
}


void PrimalDualSVM::multiply(double scalar) {
  if (scalar != 1.0) {
    w->mult(scalar);
    w2->mult(scalar); 
    if (alphas->getSize() > 0) {
      alphas->mult(scalar);
    }
    if (useBias) {
      bias += scalar;
    }
  }
}


void PrimalDualSVM::add(sgpp::base::Grid& grid, sgpp::base::DataVector& x, 
                        double alpha, size_t dataDim) {
  if (svs->getNrows() < budget) {
    // SGKernel evaluation operation
    std::unique_ptr<datadriven::OperationEvalSGKernel> opSGKernel =
      op_factory::createOperationEvalSGKernel(grid);      
	  
    sgpp::base::DataVector xTrans(grid.getSize());  
    opSGKernel->phi(x, xTrans, dataDim);

    svs->appendRow(x);
    alphas->append(alpha);  
    norms->append(xTrans.dotProduct(xTrans));

    sgpp::base::DataVector xTransCopy(xTrans); 
    xTrans.mult(alpha);
    w->add(xTrans);    
    xTransCopy.mult(std::abs(alpha));
    w2->add(xTransCopy);

    if (useBias) {
      bias += alpha;
    }
  }
}

  
/*void PrimalDualSVM::remove(size_t idx) {
  // ToDo:
}*/


}  // namespace datadriven
}  // namespace sgpp

