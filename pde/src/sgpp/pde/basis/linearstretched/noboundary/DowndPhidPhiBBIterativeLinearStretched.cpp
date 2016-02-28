// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/noboundary/DowndPhidPhiBBIterativeLinearStretched.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

DowndPhidPhiBBIterativeLinearStretched::DowndPhidPhiBBIterativeLinearStretched(
    SGPP::base::GridStorage* storage)
    : storage(storage) {}

DowndPhidPhiBBIterativeLinearStretched::~DowndPhidPhiBBIterativeLinearStretched() {}

void DowndPhidPhiBBIterativeLinearStretched::operator()(SGPP::base::DataVector& alpha,
                                                        SGPP::base::DataVector& result,
                                                        size_t dim) {
  SGPP::base::Stretching* stretching = this->storage->getStretching();
  //  float_t q = stretching->getIntervalWidth(dim);
  //
  //  float_t Qqout = 1.0/q;

  /*  // init the coefficients of the ansatz functions with boundary
    result.setAll(0.0);

    if (q != 1.0)
    {
      // traverse all basis function by sequence number
      for(size_t i = 0; i < storage->getSize(); i++)
      {
        SGPP::base::GridStorage::index_type::level_type level;
        SGPP::base::GridStorage::index_type::index_type index;
        (*storage)[i]->get(dim, level, index);
        //only affects the diagonal of the stiffness matrix
        result[i] = alpha[i]*(Qqout*pow(2.0, static_cast<int>(level+1)));
      }
    }
    else
    {
      // traverse all basis function by sequence number
      for(size_t i = 0; i < storage->getSize(); i++)
      {
        SGPP::base::GridStorage::index_type::level_type level;
        SGPP::base::GridStorage::index_type::index_type index;
        (*storage)[i]->get(dim, level, index);
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
    SGPP::base::GridStorage::index_type::level_type level;
    SGPP::base::GridStorage::index_type::index_type index;
    (*storage)[i]->get(dim, level, index);
    float_t posl = 0, posr = 0, posc = 0;
    stretching->getAdjacentPositions(static_cast<int>(level), static_cast<int>(index), dim, posc,
                                     posl, posr);
    float_t baseLength = posr - posl;
    float_t leftLength = posc - posl;
    float_t rightLength = posr - posc;
    // only affects the diagonal of the stiffness matrix
    result[i] = alpha[i] * baseLength / (leftLength * rightLength);
  }
}
}  // namespace pde
}  // namespace SGPP
