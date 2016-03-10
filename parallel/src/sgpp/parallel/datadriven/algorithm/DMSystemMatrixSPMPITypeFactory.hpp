// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP
#define DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBaseSP.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

class DMSystemMatrixSPMPITypeFactory {
 private:
  template <typename Kernel>
  static sgpp::datadriven::DMSystemMatrixBaseSP* createDMSystemMatrixMPITypeSP(
      sgpp::base::Grid& grid, sgpp::base::DataMatrixSP& trainDataset, float lambda,
      VectorizationType vecType, MPIType mpiType);

 public:
  static sgpp::datadriven::DMSystemMatrixBaseSP* getDMSystemMatrixSP(
      sgpp::base::Grid& grid, sgpp::base::DataMatrixSP& trainDataset, float lambda,
      VectorizationType vecType, MPIType mpiType);
};
}  // namespace parallel
}  // namespace sgpp

#endif  // DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP
