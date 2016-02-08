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


namespace SGPP {
namespace parallel {

class DMSystemMatrixSPMPITypeFactory {
 private:
  template<typename Kernel>
  static SGPP::datadriven::DMSystemMatrixBaseSP* createDMSystemMatrixMPITypeSP(
    SGPP::base::Grid& grid, SGPP::base::DataMatrixSP& trainDataset, float lambda,
    VectorizationType vecType, MPIType mpiType);

 public:
  static SGPP::datadriven::DMSystemMatrixBaseSP* getDMSystemMatrixSP(
    SGPP::base::Grid& grid, SGPP::base::DataMatrixSP& trainDataset, float lambda,
    VectorizationType vecType, MPIType mpiType);

};

}
}

#endif // DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP