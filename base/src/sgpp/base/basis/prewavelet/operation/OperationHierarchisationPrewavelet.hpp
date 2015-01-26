/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)


#ifndef OPERATIONHIERARCHISATIONPREWAVELET_HPP
#define OPERATIONHIERARCHISATIONPREWAVELET_HPP

#include <sgpp/base/operation/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid with prewavelets and no boundary. Please note, that there is
     * no efficient way to directly calculate the hierarchical surpluses for the prewavelet base.
     * But there is a fast and efficient way to transform hierarchical surpluses of the normal
     * linear ansatzfunctions into the prewavelet base and vice versa. Thus, we use the normal
     * hierarchisation of the linear basis and afterwards transform the resulting Vector into
     * prewavelets (see ConvertLinearToPrewavelet.hpp). For the Dehierarchisation, this process is
     * reversed (see ConvertPrewaveletToLinear.hpp)
     */
    class OperationHierarchisationPrewavelet: public OperationHierarchisation {
      public:
        /**
         * Constructor of OperationHierarchisationPrewavelet
         *
           * An adaptive grid with prewavelet ansatz functions requires for operations
         * using the up-down algorithm shadow points. These shadow points a needed just
         * for data transport, thus they do not have an influence on the final function.
         * Please refer to SGPP::pde::UpDownOneOpDimWithShadow for more information.
           *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param shadowStorage shadow points (see detailed description)
         */
        OperationHierarchisationPrewavelet(GridStorage* storage, GridStorage* shadowStorage) :
          storage(storage), shadowStorage(shadowStorage) {
        }

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationPrewavelet() {
        }

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
        GridStorage* shadowStorage;

        void expandGrid();

        void shrinkGrid();
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONPREWAVELET_HPP */
