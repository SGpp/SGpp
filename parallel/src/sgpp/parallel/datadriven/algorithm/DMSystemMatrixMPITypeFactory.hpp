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

namespace SGPP {
namespace parallel {

class DMSystemMatrixMPITypeFactory {
 private:
  template <typename Kernel>
  static SGPP::datadriven::DMSystemMatrixBase* createDMSystemMatrixMPIType(
      SGPP::base::Grid& grid, SGPP::base::DataMatrix& trainDataset, double lambda,
      VectorizationType vecType, MPIType mpiType);

 public:
  static SGPP::datadriven::DMSystemMatrixBase* getDMSystemMatrix(
      SGPP::base::Grid& grid, SGPP::base::DataMatrix& trainDataset, double lambda,
      VectorizationType vecType, MPIType mpiType);
};

}  // namespace parallel
}  // namespace SGPP

#endif  // DMSYSTEMMATRIXMPITYPEFACTORY_HPP
