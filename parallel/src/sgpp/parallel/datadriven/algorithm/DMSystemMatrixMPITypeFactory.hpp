/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXMPITYPEFACTORY_HPP
#define DMSYSTEMMATRIXMPITYPEFACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

namespace sg {
  namespace parallel {

    class DMSystemMatrixMPITypeFactory {
      private:
        template<typename Kernel>
        static sg::datadriven::DMSystemMatrixBase* createDMSystemMatrixMPIType(sg::base::Grid& grid, sg::base::DataMatrix& trainDataset, double lambda, VectorizationType vecType, MPIType mpiType);

      public:
        static sg::datadriven::DMSystemMatrixBase* getDMSystemMatrix(sg::base::Grid& grid, sg::base::DataMatrix& trainDataset, double lambda, VectorizationType vecType, MPIType mpiType);
    };

  }
}

#endif // DMSYSTEMMATRIXMPITYPEFACTORY_HPP
