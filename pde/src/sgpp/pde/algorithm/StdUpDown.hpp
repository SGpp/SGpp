// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef STDUPDOWN_HPP
#define STDUPDOWN_HPP

#include <vector>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#ifndef TASKS_PARALLEL_UPDOWN
#define TASKS_PARALLEL_UPDOWN 4
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implements a standard Up/Down Schema without any operation dim.
     *
     * @version $HEAD$
     */
    class StdUpDown: public SGPP::base::OperationMatrix {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        StdUpDown(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~StdUpDown();


        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        /**
         * this functions provides the same functionality as the normal mult routine.
         * However, it doesn't set up an OpenMP task initialization as the mult routine.
         * This method has to be called within a OpenMP task parallelized region.
         *
         * Using this function is useful in following case: Assuming the solver of a certain
         * requires several operators in the space discretization (e.g. Black Scholes Equations)
         * this method can be used to parallelize their calculation which might results results
         * in a better parallel efficiency on systems with 4 or more cores hence fewer barriers
         * are needed.
         *
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        void multParallelBuildingBlock(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);


      protected:
        typedef SGPP::base::GridStorage::grid_iterator grid_iterator;

        /// Pointer to the grid's storage object
        SGPP::base::GridStorage* storage;
        /// algorithmic dimensions, operator is applied in this dimensions
        const std::vector<size_t> algoDims;
        /// number of algorithmic dimensions
        const size_t numAlgoDims_;
        /// max number of parallel stages (dimension recursive calls)
        static const size_t maxParallelDims_ = TASKS_PARALLEL_UPDOWN;

        /**
         * Recursive procedure for updown
         *
         * @param dim the current dimension
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        void updown(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * 1D up Operation
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down Operation
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) = 0;
    };

  }
}

#endif /* STDUPDOWN_HPP */