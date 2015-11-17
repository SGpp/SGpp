// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONSTENCILHIERARCHISATIONLINEAR_HPP
#define OPERATIONSTENCILHIERARCHISATIONLINEAR_HPP

#include <sgpp/base/operation/hash/OperationStencilHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, linear grid without boundaries
     */
    class OperationStencilHierarchisationLinear : public OperationStencilHierarchisation {
      public:
        /**
         * Constructor of OperationStencilHierarchisationLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationStencilHierarchisationLinear(GridStorage* storage) : storage(storage),
          surplusStencil(0), neighborStencil(0), weightStencil(0) {}

        /**
         * Destructor
         */
        virtual ~OperationStencilHierarchisationLinear() override {}

        virtual void doHierarchisation(DataVector& node_values) override;
        virtual void doDehierarchisation(DataVector& alpha) override;


        virtual const IndexStencil&
        getSurplusStencil() const override {
          return surplusStencil;
        };

        virtual const IndexStencil&
        getNeighborStencil() const override {
          return neighborStencil;
        };

        virtual const WeightStencil&
        getWeightStencil() const override {
          return weightStencil;
        };

        virtual size_t
        getStencilSize() const override {
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

#endif /* OPERATIONSTENCILHIERARCHISATION_HPP */
