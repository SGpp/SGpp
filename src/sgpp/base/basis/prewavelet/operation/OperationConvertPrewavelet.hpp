/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONCONVERTPREWAVELET_HPP
#define OPERATIONCONVERTPREWAVELET_HPP

#include "base/operation/OperationConvert.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     *
     */
    class OperationConvertPrewavelet : public OperationConvert {
      public:
        /**
         * Constructor of OperationHierarchisationPrewavelet
         *
           * An adaptive grid with prewavelet ansatz functions requires for operations
         * using the up-down algorithm shadow points. These shadow points a needed just
         * for data transport, thus they do not have an influence on the final function.
         * Please refer to sg::pde::UpDownOneOpDimWithShadow for more information.
         *
           * @param storage Pointer to the grid's gridstorage obejct
         * @param shadowstorage shadow points (see detailed description)
         */
        OperationConvertPrewavelet(GridStorage* storage, GridStorage* shadowstorage) :
          storage(storage), shadowstorage(shadowstorage) {
        }

        /**
         * Destructor
         */
        virtual ~OperationConvertPrewavelet() {
        }

        virtual void doConvertToLinear(DataVector& alpha);
        virtual void doConvertFromLinear(DataVector& alpha);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
        GridStorage* shadowstorage;
    };

  }
}

#endif /* OPERATIONCONVERTPREWAVELET_HPP */
