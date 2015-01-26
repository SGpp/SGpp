// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAR_HPP
#define OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAR_HPP

#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>


#include <sgpp/parallel/tools/TypesParallel.hpp>

#if defined(__SSE4_2__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /**
    * Implementation for linear functions of LTwoDotProduct Operation, linear grids without boundaries
     *
     * @version $HEAD$
     */
    class OperationLTwoDotProductVectorizedLinear: public SGPP::base::OperationMatrix {
      private:

        SGPP::base::GridStorage* storage;
        SGPP::base::DataMatrix* level_;
        SGPP::base::DataMatrix* level_int_;
        SGPP::base::DataMatrix* index_;
        SGPP::base::DataVector* lcl_q_;
        SGPP::base::DataVector* alpha_padded_;
        SGPP::base::DataVector* constants_;

        SGPP::base::DataVector** gradient_temp;
        SGPP::base::DataVector** l2dot_temp;

#if defined(STORE_MATRIX)
        SGPP::base::DataMatrix* operation_result_matrix_;
        bool operation_result_generated_;
#endif

        int process_count;
        int process_index;

        std::vector<int> all_i_start;
        std::vector<int> all_i_size;

        std::vector<int> send_start;
        std::vector<int> send_size;

        std::vector<int> recv_start;
        std::vector<int> recv_size;


        void init_constants();
        void init_grid_storage();

        double l2dot(size_t i, size_t j, size_t dim);
        double all_time;
        double all_iterations;
        SGPP::base::SGppStopwatch stopWatch;

      public:
        /**
         * Construtor of OperationLTwoDotProductVectorizedLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLTwoDotProductVectorizedLinear(SGPP::base::GridStorage* storage);

        /**
         * Construtor of OperationLTwoDotProductVectorizedLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda Vector which contains pre-factors for every dimension of the operator
         */
        OperationLTwoDotProductVectorizedLinear(SGPP::base::GridStorage* storage, SGPP::base::DataVector& lambda);

        /**
         * Destructor
         */
        virtual ~OperationLTwoDotProductVectorizedLinear();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        virtual void reset();
    };

  }

}

#endif /* OPERATIONLAPLACEVECTORIZEDLINEAR_HPP */