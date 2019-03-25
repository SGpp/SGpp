// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixDRE.hpp>
// #include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

DMSystemMatrixDRE::DMSystemMatrixDRE(sgpp::base::Grid& grid,
                                     sgpp::base::DataMatrix& lhsData,
                                     sgpp::base::DataMatrix& rhsData,
                                     std::shared_ptr<base::OperationMatrix> C,
                                     double lambdaRegression)
    : DMSystemMatrixBase(lhsData, lambdaRegression),
      grid(grid),
      rhs_dataset_(rhsData),
      C(std::move(C)) {
  // this->B = sgpp::op_factory::createOperationMultiEval(grid);
  this->B.reset(
      sgpp::op_factory::createOperationMultipleEval(grid, this->dataset_));
}

DMSystemMatrixDRE::~DMSystemMatrixDRE() {}

void DMSystemMatrixDRE::mult(sgpp::base::DataVector& alpha,
                             sgpp::base::DataVector& result) {
  // temp variable may be resized during computation
  sgpp::base::DataVector temp(this->dataset_.getNrows());
  size_t Mp = this->dataset_.getNrows();

  // Operation B
  //      this->B->mult((*this->dataset_), alpha, temp);
  //      this->B->multTranspose((*this->dataset_), temp, result);

  // this->B->mult(alpha, temp);

  std::unique_ptr<base::OperationMultipleEval> op(
      sgpp::op_factory::createOperationMultipleEval(grid, this->dataset_));
  // sgpp::base::DataVector temp2(this->dataset_.getNrows());
  // op->mult(*(this->dataset_), alpha, temp2);
  op->mult(alpha, temp);

  /*if (temp.getSize() != temp2.getSize()) {
   std::cout << "error: sizes don't match" << std::endl;
   std::cout << "ref: " << temp.getSize() << " mine: " << temp2.getSize() <<
   std::endl;
   }*/

  /*for (size_t i = 0; i < temp.getSize(); i++) {
   if (temp[i] != temp2[i]) {
   if (fabs(temp[i] - temp2[i]) > 0.0000001) {
   std::cout << "ref: " << temp[i] << " mine: " << temp2[i] << " i: " << i <<
   std::endl;
   }
   }
   }*/

  //      this->B->multTranspose(temp, result);
  //      sgpp::base::DataVector result2(result);
  //      op->multTranspose(*(this->dataset_), temp, result2);
  //      if (result.getSize() != result2.getSize()) {
  //        std::cout << "ref: " << result.getSize() << " mine: " <<
  //        result2.getSize() << std::endl;
  //      }
  //      for (size_t i = 0; i < result.getSize(); i++) {
  //        if (result[i] != result2[i]) {
  //          if (fabs(result[i] - result2[i]) > 0.0000001) {
  //            std::cout << "ref: " << result[i] << " mine: "  << result2[i] <<
  //            " i: " << i <<
  //              std::endl;
  //          }
  //        }
  //      }
  op->multTranspose(temp, result);

  sgpp::base::DataVector temptwo(alpha.getSize());
  this->C->mult(alpha, temptwo);
  result.axpy(static_cast<double>(Mp) * this->lambda_, temptwo);
}

void DMSystemMatrixDRE::generateb(sgpp::base::DataVector& b) {
  // this->B->multTranspose((*this->dataset_), classes, b);
  // this->B->multTranspose(classes, b);

  // rhs dataset is used here, weighted by the ratio of dataset sizes
  std::unique_ptr<base::OperationMultipleEval> op(
      sgpp::op_factory::createOperationMultipleEval(grid, this->rhs_dataset_));
  double Mp = static_cast<double>(this->dataset_.getNrows());
  double Mq = static_cast<double>(this->rhs_dataset_.getNrows());
  sgpp::base::DataVector y(b.getSize());
  y.setAll(Mp / Mq);
  op->multTranspose(y, b);
}

}  // namespace datadriven
}  // namespace sgpp
