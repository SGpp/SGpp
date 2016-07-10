// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/noboundary/DowndPhidPhiBBIterativeLinearStretched.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

DowndPhidPhiBBIterativeLinearStretched::DowndPhidPhiBBIterativeLinearStretched(
    sgpp::base::GridStorage* storage)
    : storage(storage) {}

DowndPhidPhiBBIterativeLinearStretched::~DowndPhidPhiBBIterativeLinearStretched() {}

void DowndPhidPhiBBIterativeLinearStretched::operator()(sgpp::base::DataVector& alpha,
                                                        sgpp::base::DataVector& result,
                                                        size_t dim) {
  sgpp::base::Stretching* stretching = this->storage->getStretching();
  //  double q = stretching->getIntervalWidth(dim);
  //
  //  double Qqout = 1.0/q;

  /*  // init the coefficients of the ansatz functions with boundary
    result.setAll(0.0);

    if (q != 1.0)
    {
      // traverse all basis function by sequence number
      for(size_t i = 0; i < storage->getSize(); i++)
      {
        sgpp::base::level_t level;
        sgpp::base::index_t index;
        (*storage)[i].get(dim, level, index);
        //only affects the diagonal of the stiffness matrix
        result[i] = alpha[i]*(Qqout*pow(2.0, static_cast<int>(level+1)));
      }
    }
    else
    {
      // traverse all basis function by sequence number
      for(size_t i = 0; i < storage->getSize(); i++)
      {
        sgpp::base::level_t level;
        sgpp::base::index_t index;
        (*storage)[i].get(dim, level, index);
        //only affects the diagonal of the stiffness matrix
        result[i] = alpha[i]*pow(2.0, static_cast<int>(level+1));
      }
    }*/
  //  else{
  //    std::cout<<"else called"<<std::endl;
  //  }
  result.setAll(0.0);

  // traverse all basis function by sequence number
  for (size_t i = 0; i < storage->getSize(); i++) {
    sgpp::base::level_t level;
    sgpp::base::index_t index;
    (*storage)[i].get(dim, level, index);
    double posl = 0, posr = 0, posc = 0;
    stretching->getAdjacentPositions(static_cast<int>(level), static_cast<int>(index), dim, posc,
                                     posl, posr);
    double baseLength = posr - posl;
    double leftLength = posc - posl;
    double rightLength = posr - posc;
    // only affects the diagonal of the stiffness matrix
    result[i] = alpha[i] * baseLength / (leftLength * rightLength);
  }
}
}  // namespace pde
}  // namespace sgpp
