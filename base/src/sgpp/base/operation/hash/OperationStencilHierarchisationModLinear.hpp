// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONSTENCILHIERARCHISATIONMODLINEAR_HPP
#define OPERATIONSTENCILHIERARCHISATIONMODLINEAR_HPP

#include <sgpp/base/operation/hash/OperationStencilHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, linear grid with modified basis functions
     */
    class OperationStencilHierarchisationModLinear : public OperationStencilHierarchisation {
      public:
        /**
         * Constructor of OperationStencilHierarchisationModLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationStencilHierarchisationModLinear(GridStorage* storage) : storage(storage),
          surplusStencil(0), neighborStencil(0), weightStencil(0) {}

        /**
         * Destructor
         */
        virtual ~OperationStencilHierarchisationModLinear() {}

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);


        virtual const IndexStencil&
        getSurplusStencil() const {
          return surplusStencil;
        };

        virtual const IndexStencil&
        getNeighborStencil() const {
          return neighborStencil;
        };

        virtual const WeightStencil&
        getWeightStencil() const  {
          return weightStencil;
        };

        virtual size_t
        getStencilSize() const {
          return surplusStencil.size();
        };

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;

        /// Index array with surplus indices
        IndexStencil surplusStencil;

        /// Index array with neighboring surplus indices
        IndexStencil neighborStencil;

        /// Index array with surplus indices
        WeightStencil weightStencil;

    };

  }
}

#endif /* OPERATIONSTENCILHIERARCHISATIONMODLINEAR_HPP */