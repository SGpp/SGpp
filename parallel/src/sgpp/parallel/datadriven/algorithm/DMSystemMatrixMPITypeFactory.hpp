// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXMPITYPEFACTORY_HPP
#define DMSYSTEMMATRIXMPITYPEFACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

class DMSystemMatrixMPITypeFactory {
 private:
  template <typename Kernel>
  static sgpp::datadriven::DMSystemMatrixBase* createDMSystemMatrixMPIType(
      sgpp::base::Grid& grid, sgpp::base::DataMatrix& trainDataset, double lambda,
      VectorizationType vecType, MPIType mpiType);

 public:
  static sgpp::datadriven::DMSystemMatrixBase* getDMSystemMatrix(
      sgpp::base::Grid& grid, sgpp::base::DataMatrix& trainDataset, double lambda,
      VectorizationType vecType, MPIType mpiType);
};

}  // namespace parallel
}  // namespace sgpp

#endif  // DMSYSTEMMATRIXMPITYPEFACTORY_HPP
