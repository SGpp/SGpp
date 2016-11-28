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

// -------------------- constructors and destructors --------------------
PrimalDualSVM::PrimalDualSVM(size_t dim, size_t dataDim, int B, bool useBias) // dim = size of grid, dataDim = dimension of input data
    : svs(new sgpp::base::DataMatrix(0,dataDim)),       
      alphas(new sgpp::base::DataVector(0)),
      norms(new sgpp::base::DataVector(0)), 
      w(new sgpp::base::DataVector(dim, 0.0)),    // dim of transformed feature space (grid size), w needs to be resized if grid gets resized !!!
      w2(new sgpp::base::DataVector(dim, 0.0)), 
      b(0.0),
      B(B),
      useBias(useBias) {

  svs->setInc(B);
}

PrimalDualSVM::~PrimalDualSVM() {}
// -----------------------------------------------------------------------


double PrimalDualSVM::predictRaw(sgpp::base::Grid& grid, sgpp::base::DataVector& x, size_t dataDim, bool trans) {
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
    res += b;
  }
  return res;
}


int PrimalDualSVM::predict(sgpp::base::Grid& grid, sgpp::base::DataVector& x, size_t dataDim) {
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
      b += scalar;
    }
  }
}

void PrimalDualSVM::add(sgpp::base::Grid& grid, sgpp::base::DataVector& x, double alpha, size_t dataDim) {
  //SGKernel evaluation operation
  std::unique_ptr<datadriven::OperationEvalSGKernel> opSGKernel =
    op_factory::createOperationEvalSGKernel(grid);      
	  
  sgpp::base::DataVector xTrans(grid.getSize());  
  opSGKernel->phi(x, xTrans, dataDim);

  svs->appendRow(x);          // svs = dataMatrix
  alphas->append(alpha);  
  norms->append(xTrans.dotProduct(xTrans));

  sgpp::base::DataVector xTransCopy(xTrans); 
  xTrans.mult(alpha);
  w->add(xTrans);    
  xTransCopy.mult(std::abs(alpha));
  w2->add(xTransCopy);

  if (useBias) {
    b += alpha;
  }
}
  
/*void PrimalDualSVM::remove(size_t idx) {
//ToDo:
}*/

}  // namespace datadriven
}  // namespace sgpp
