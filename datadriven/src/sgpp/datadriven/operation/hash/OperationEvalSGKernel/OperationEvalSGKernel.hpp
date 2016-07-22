// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OperationEvalSGKernel_HPP_
#define OperationEvalSGKernel_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class OperationEvalSGKernel {
 protected:
  sgpp::base::Grid& grid;
  
 public:
  OperationEvalSGKernel(base::Grid& grid);
  
  virtual ~OperationEvalSGKernel();
  
  virtual void phi(sgpp::base::DataVector& x, sgpp::base::DataVector& xTrans, size_t dataDim);
  
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* OperationEvalSGKernel_HPP_ */
