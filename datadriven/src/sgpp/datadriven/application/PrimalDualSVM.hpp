// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PRIMALDUALSVM_HPP_
#define PRIMALDUALSVM_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

class PrimalDualSVM {
 protected:

  double b; 
  
  int B;
  
  bool useBias;
 
 public:

  std::shared_ptr<base::DataMatrix> svs; 
  
  std::shared_ptr<base::DataVector> alphas; 
  
  std::shared_ptr<base::DataVector> norms; 
  
  std::shared_ptr<base::DataVector> w; 
  
  std::shared_ptr<base::DataVector> w2;

  PrimalDualSVM(size_t dim, size_t inputDim, int B, bool useBias);
  
  virtual ~PrimalDualSVM();
  
  virtual double predictRaw(sgpp::base::Grid& grid, sgpp::base::DataVector& x, size_t dataDim, bool trans = false);
  
  virtual int predict(sgpp::base::Grid& grid, sgpp::base::DataVector& x, size_t dataDim);
  
  virtual void add(sgpp::base::Grid& grid, sgpp::base::DataVector& x, double alpha, size_t dataDim); 
  
  virtual void multiply(double scalar);
  
  virtual void remove(size_t idx); 
  
  //virtual void maintainSVS();
  
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* PRIMALDUALSVM_HPP_ */
