/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

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
        static SGPP::datadriven::DMSystemMatrixBaseSP* createDMSystemMatrixMPITypeSP(SGPP::base::Grid& grid, SGPP::base::DataMatrixSP& trainDataset, float lambda, VectorizationType vecType, MPIType mpiType);

      public:
        static SGPP::datadriven::DMSystemMatrixBaseSP* getDMSystemMatrixSP(SGPP::base::Grid& grid, SGPP::base::DataMatrixSP& trainDataset, float lambda, VectorizationType vecType, MPIType mpiType);

    };

  }
}

#endif // DMSYSTEMMATRIXSPMPITYPEFACTORY_HPP
