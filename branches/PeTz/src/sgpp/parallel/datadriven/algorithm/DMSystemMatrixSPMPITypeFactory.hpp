/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP
#define DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP

#include "base/grid/Grid.hpp"
#include "parallel/tools/TypesParallel.hpp"
#include "datadriven/algorithm/DMSystemMatrixBaseSP.hpp"

namespace sg {
  namespace parallel {

    class DMSystemMatrixSPMPITypeFactory {
      private:
        template<typename Kernel>
        static sg::datadriven::DMSystemMatrixBaseSP* createDMSystemMatrixMPITypeSP(sg::base::Grid& grid, sg::base::DataMatrixSP& trainDataset, float lambda, VectorizationType vecType, MPIType mpiType);

      public:
        static sg::datadriven::DMSystemMatrixBaseSP* getDMSystemMatrixSP(sg::base::Grid& grid, sg::base::DataMatrixSP& trainDataset, float lambda, VectorizationType vecType, MPIType mpiType);

    };

  }
}

#endif // DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP
