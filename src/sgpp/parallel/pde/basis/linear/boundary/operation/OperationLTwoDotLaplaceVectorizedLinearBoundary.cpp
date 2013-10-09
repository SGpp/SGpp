/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)


#ifdef USE_MPI
#include <mpi.h>
#endif

#include <omp.h>

#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/exception/operation_exception.hpp"

#include "parallel/pde/basis/linear/boundary/operation/OperationLTwoDotLaplaceVectorizedLinearBoundary.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/tools/TypesParallel.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "base/tools/SGppStopwatch.hpp"

#include <cmath>
#include <assert.h>


#include <cstring>
#include <iostream>
#include <vector>

#if defined(__SSE4_2__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

#if defined(__SSE4_2__) && !defined(__AVX__)
#define VECTOR_SIZE 2
#else
#define VECTOR_SIZE 4
#endif


#define REG_BCOUNT 2

#define BLOCK_LENGTH (REG_BCOUNT * VECTOR_SIZE)
#define PAGE_TRUE_CAPACITY 512

#define ROUND_DOWN(X, Y) (((X) / (Y)) * (Y))
#define ROUND_UP(X, Y) (std::ceil((1.0 * X) / (Y)) * (Y))


//#define WITH_GLFOPS_CALCULATION


namespace sg {
    namespace parallel {
			
        OperationLTwoDotLaplaceVectorizedLinearBoundary::OperationLTwoDotLaplaceVectorizedLinearBoundary(sg::base::GridStorage* storage) 
			: 
			storage(storage),
             level_(NULL),
             level_int_(NULL),
             index_(NULL),
             lcl_q_(NULL),
             lcl_q_inv_(NULL),
             lambda_(NULL),
             constants_(NULL),
			 alpha_padded_(NULL),
			 result_boundary_filtered_(NULL),
			 level_boundary_filtered_(NULL),
             level_int_boundary_filtered_(NULL),
             index_boundary_filtered_(NULL),
			gradient_temp(NULL),
			l2dot_temp(NULL)
#if defined(STORE_MATRIX)
		,
		operation_result_matrix_(NULL)
#endif
		{
std::cout<<"IN CONSTRUCTOR: OperationLTwoDotLaplaceVectorizedLinearBoundary" << std::endl;
            init_constants();
			init_grid_storage();
			
			this->lambda_ = new sg::base::DataVector(storage->dim());
            this->lambda_->setAll(1.0);
        }
        
        OperationLTwoDotLaplaceVectorizedLinearBoundary::OperationLTwoDotLaplaceVectorizedLinearBoundary(sg::base::GridStorage* storage, sg::base::DataVector& lambda)
			: 
			storage(storage),
             level_(NULL),
             level_int_(NULL),
             index_(NULL),
             lcl_q_(NULL),
             lcl_q_inv_(NULL),
             lambda_(NULL),
             constants_(NULL),
			 alpha_padded_(NULL),
			 result_boundary_filtered_(NULL),
			 level_boundary_filtered_(NULL),
             level_int_boundary_filtered_(NULL),
             index_boundary_filtered_(NULL),
			gradient_temp(NULL),
			l2dot_temp(NULL)
#if defined(STORE_MATRIX)
		,
		operation_result_matrix_(NULL)
#endif
		{
std::cout<<"IN CONSTRUCTOR: OperationLTwoDotLaplaceVectorizedLinearBoundary" << std::endl;
            init_constants();
			init_grid_storage();
			
            this->lambda_ = new sg::base::DataVector(lambda);
        }
        
        OperationLTwoDotLaplaceVectorizedLinearBoundary::~OperationLTwoDotLaplaceVectorizedLinearBoundary() {
            delete this->level_;
            delete this->level_int_;
            delete this->index_;
            delete lcl_q_;
            delete lcl_q_inv_;
			delete constants_;
			
            delete this->alpha_padded_;
            
            if (this->lambda_ != NULL)
                delete this->lambda_;

			 delete result_boundary_filtered_;
			 delete level_boundary_filtered_;
             delete level_int_boundary_filtered_;
             delete index_boundary_filtered_;
			 delete gradient_temp;
			 delete l2dot_temp;
#if defined(STORE_MATRIX)
			 delete operation_result_matrix_;
#endif
        }
		
		void OperationLTwoDotLaplaceVectorizedLinearBoundary::reset(){
			init_grid_storage();
		}        
		
		void OperationLTwoDotLaplaceVectorizedLinearBoundary::init_constants(){		
            #if defined(__SSE4_2__)
            this->constants_ = new sg::base::DataVector(0);
            
            this->constants_->append(0);
            this->constants_->append(0.5);
            this->constants_->append(2.0 / 3.0);
            this->constants_->append(1.0);
            this->constants_->append(2.0);
            
            unsigned long long int abs_mask = 0x7fffffffffffffff;
            this->constants_->append(*(reinterpret_cast<double*>(&abs_mask)));
            #endif
			
			process_index = 0;
			process_count = 1;
#ifdef USE_MPI 
            MPI_Comm_rank(MPI_COMM_WORLD, &process_index);
            MPI_Comm_size(MPI_COMM_WORLD, &process_count);
#endif	
		}
		
        void OperationLTwoDotLaplaceVectorizedLinearBoundary::init_grid_storage() {
			if(this->level_)
				delete this->level_;
            this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
            
			if(this->level_int_)
				delete this->level_int_;
			this->level_int_ = new sg::base::DataMatrix(storage->size(), storage->dim());
			
			if(this->index_)
				delete this->index_;
			this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());
			
			if(this->lcl_q_)
				delete this->lcl_q_;
            this->lcl_q_ = new sg::base::DataVector(this->storage->dim());
            double* lcl_q_ptr_ = lcl_q_->getPointer();
			
			if(this->lcl_q_inv_)
				delete this->lcl_q_inv_;
            this->lcl_q_inv_ = new sg::base::DataVector(this->storage->dim());
            double* lcl_q_inv_ptr_ = lcl_q_inv_->getPointer();
            
            storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
            storage->getLevelForIntegral(*(this->level_int_));
            
            
            for(size_t i = 0; i < this->storage->size(); i++)
            {
                bool i_boundary = false;
                
                for (size_t d = 0; d < this->storage->dim(); d++) {
                    i_boundary = i_boundary || (level_->get(i, d) == 1);
                }
                
                if (i_boundary == false)
                    i_boundary_filtered.push_back(i);
            }
			
			
			if(this->result_boundary_filtered_)
				delete this->result_boundary_filtered_;
            this->result_boundary_filtered_ = new sg::base::DataVector((size_t) i_boundary_filtered.size());
			
			if(level_boundary_filtered_)
				delete level_boundary_filtered_;
            level_boundary_filtered_ = new sg::base::DataMatrix((size_t) i_boundary_filtered.size(), level_->getNcols());
            
			if(level_int_boundary_filtered_)
				delete level_int_boundary_filtered_;
			level_int_boundary_filtered_ = new sg::base::DataMatrix((size_t) i_boundary_filtered.size(), level_int_->getNcols());
            
			if(index_boundary_filtered_)
				delete index_boundary_filtered_;
			index_boundary_filtered_ = new sg::base::DataMatrix((size_t) i_boundary_filtered.size(), index_->getNcols());
			
			//used for copying rows
            sg::base::DataVector level_copy_vector(this->level_->getNcols());
            sg::base::DataVector level_int_copy_vector(this->level_int_->getNcols());
            sg::base::DataVector index_copy_vector(this->index_->getNcols());

			for(size_t i = 0; i < i_boundary_filtered.size(); ++i)
            {
                level_->getRow(i_boundary_filtered[i], level_copy_vector);
                level_boundary_filtered_->setRow(i, level_copy_vector);
                
                level_int_->getRow(i_boundary_filtered[i], level_int_copy_vector);
                level_int_boundary_filtered_->setRow(i, level_int_copy_vector);
                
                index_->getRow(i_boundary_filtered[i], index_copy_vector);
                index_boundary_filtered_->setRow(i, index_copy_vector);
			}
			
			storage->getLevelIndexArraysForEvalTLBOptimized(*(this->level_), *(this->index_), sg::parallel::X86SIMD, BLOCK_LENGTH);
            storage->getLevelForIntegralTLBOptimized(*(this->level_int_), sg::parallel::X86SIMD, BLOCK_LENGTH);
			
			std::size_t padded_size = this->level_->getNcols();
			if(alpha_padded_)
				delete alpha_padded_;
            alpha_padded_ = new sg::base::DataVector(padded_size);
            alpha_padded_->setAll(0.0);
			
            size_t single_process_portion = (this->i_boundary_filtered.size() / process_count) + 1;

			all_i_start.clear();
			all_i_size.clear();
			send_start.clear();
			send_size.clear();
			recv_start.clear();
			recv_size.clear();
			
            for(int i = 0; i < process_count; ++i)
            {
				int process_start = (int) (i * single_process_portion);
                all_i_start.push_back(process_start);
                int process_portion = (i == process_count - 1)? (int) (this->i_boundary_filtered.size() - i * single_process_portion)
																: (int) (single_process_portion);
                
				process_portion = std::min<int>(process_portion, (int) (this->i_boundary_filtered.size() - i * single_process_portion));
				process_portion = std::max<int>(process_portion, 0);
				
				all_i_size.push_back(process_portion);
            }
            
            for(int i = 0; i < process_count; ++i)
            {
                send_start.push_back((int) (single_process_portion * process_index));
                int process_send_size = (process_index == process_count -1)?  (int) (this->i_boundary_filtered.size() - single_process_portion * process_index)
																			: (int) (single_process_portion);
																			
                process_send_size = std::min<int>(process_send_size, (int) (this->i_boundary_filtered.size() - single_process_portion * process_index));
				process_send_size = std::max<int>(process_send_size, 0);
				//process_send_size = (i == process_index)? 0: process_send_size;
				
                send_size.push_back(process_send_size);
                
                recv_start.push_back((int) (single_process_portion * i));
                int process_recv_size = (i == process_count -1)?  (int) (this->i_boundary_filtered.size() - single_process_portion * i)
																: (int) (single_process_portion);

                process_recv_size = std::min<int>(process_recv_size, (int) (this->i_boundary_filtered.size() - single_process_portion * i));
				process_recv_size = std::max<int>(process_recv_size, 0);
				//process_recv_size = (i == process_index)? 0 : process_recv_size;
				
                recv_size.push_back(process_recv_size);
            }
			
#pragma omp parallel
{
				if(this->gradient_temp)
					delete this->gradient_temp[omp_get_thread_num()];
				if(this->l2dot_temp)
					delete this->l2dot_temp[omp_get_thread_num()];					
#pragma omp barrier

#pragma omp single
{
				if(this->gradient_temp)
					delete this->gradient_temp;
				if(this->l2dot_temp)
					delete this->l2dot_temp;

			this->gradient_temp = new sg::base::DataVector*[omp_get_num_threads()];
			this->l2dot_temp = new sg::base::DataVector*[omp_get_num_threads()];
}
#pragma omp barrier

			//std::cout << "OMP THREAD :" << omp_get_thread_num() << std::endl;

			gradient_temp[omp_get_thread_num()] = new sg::base::DataVector(VECTOR_SIZE * this->storage->dim() * REG_BCOUNT);
			l2dot_temp[omp_get_thread_num()] = new sg::base::DataVector(VECTOR_SIZE * this->storage->dim() * REG_BCOUNT);
}
            // fill q array
            for (size_t d = 0; d < this->storage->dim(); d++) {
                sg::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
                lcl_q_ptr_[d] = boundingBox->getIntervalWidth(d);
                lcl_q_inv_ptr_[d] = 1.0 / boundingBox->getIntervalWidth(d);
            }

			
#if defined (STORE_MATRIX)
			size_t result_matrix_rows = all_i_size[process_index];
			size_t result_matrix_cols = storage->size(); //(size_t) (ceil((double)storage->size() / BLOCK_LENGTH) * BLOCK_LENGTH);
			
			//check if matrix fits in memory 
			char* matrix_max_size_gb_str = getenv("SGPP_PDE_MATRIX_SIZE_GB");
			size_t matrix_max_size_gb = 0;
			if (matrix_max_size_gb_str == NULL || strlen(matrix_max_size_gb_str) == 0)
			{
				matrix_max_size_gb = 4;
			}
			else
			{
				std::stringstream temp_stream(matrix_max_size_gb_str);
				temp_stream >> matrix_max_size_gb;
			}
			
			size_t matrix_needed_size_bytes = result_matrix_rows * result_matrix_cols * sizeof(double);
			
			if(matrix_needed_size_bytes > matrix_max_size_gb * 1024 * 1024 * 1024)
			{
				size_t matrix_needed_size_gb = (size_t) (ceil((double)matrix_needed_size_bytes / (1024 * 1024 * 1024)));
				char exception_string[512];
				sprintf(exception_string, "OperationLaplaceVectorizedLinearBoundary::init : More memory (= %i GB) needed to store the operation matrix, Please set the SGPP_PDE_MATRIX_SIZE_GB environment variable!", (int) matrix_needed_size_gb);
				
				std::cerr << exception_string << std::endl;
		        throw new sg::base::operation_exception(exception_string);
			}

			
			if(operation_result_matrix_)
				delete operation_result_matrix_;
			
			operation_result_matrix_ = new sg::base::DataMatrix(result_matrix_rows, result_matrix_cols);
			
			//create transposed, pad, then re-transpose
			sg::parallel::VectorizationType vectorizationType = sg::parallel::X86SIMD;
			operation_result_matrix_ = new sg::base::DataMatrix(result_matrix_cols, result_matrix_rows);
			sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->operation_result_matrix_), vectorizationType);
			operation_result_matrix_->transpose();
			
			operation_result_generated_ = false;
			

#pragma omp parallel
{
		size_t padded_size = this->operation_result_matrix_->getNcols();
		size_t thr_start;
		size_t thr_end;
		sg::parallel::PartitioningTool::getOpenMPPartitionSegment(0, result_matrix_rows, &thr_start, &thr_end);
		
		for(size_t i = thr_start; i < thr_end; i++)
		{
			double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (i) * operation_result_matrix_->getNcols();
			
			for(size_t j = 0; j < padded_size; ++j)
			{
				operation_result_dest_ptr[j] = 0.0;
			}
		}
}			
#endif

        }
        
        double OperationLTwoDotLaplaceVectorizedLinearBoundary::gradient_dirichlet(size_t i, size_t j, size_t dim) {
            double grad;
            
            double i_level_grad = level_->get(i, dim);
            double i_index_grad = index_->get(i, dim);
            double j_level_grad = level_->get(j, dim);
            double j_index_grad = index_->get(j, dim);
            
            // only affects the diagonal of the stiffness matrix, on level 0 we have zero
            // due to Dirichlet boundary conditions
            bool doGrad = ((i_level_grad == j_level_grad) && (i_index_grad == j_index_grad) && (i_level_grad != 1));
            grad = i_level_grad * 2.0 * doGrad;
            
            // scale by 1/q
            return (grad * *(lcl_q_inv_->getPointer() + dim));
        }
        
        
        double OperationLTwoDotLaplaceVectorizedLinearBoundary::l2dot_dirichlet(size_t i, size_t j, size_t dim) {
            double lid = level_->get(i, dim);
            double ljd = level_->get(j, dim);
            double iid = index_->get(i, dim);
            double ijd = index_->get(j, dim);
            double in_lid = level_int_->get(i, dim);
            double in_ljd = level_int_->get(j, dim);
            
            //////////
            /// First case: lid == ljd
            /// we mask this in the end in the last line
            //////////
            
            // use textbook formular if both operands are identical
            // ansatz function on the same level but with different indecies
            // don't overlap!
            double res_one = (2.0 / 3.0) * in_lid * ((iid == ijd) && (ljd != 1));
            
            //////////
            /// Second case: lid != ljd
            /// we mask this in the end in the last line
            //////////
            
            // now we select the 1st as the "narrow" basisfunction (it has a higher level)
            // --> we know can regard 2nd function as linear function and therefore
            // apply the wellknown formular: 1/2 * (f_l + f_r) * 2^{-l}
            bool selector = (lid > ljd);
            double i1d = iid * (selector) + ijd * (!selector);
            //double l1d = lid*(selector) + ljd*(!selector);
            double in_l1d = in_lid * (selector) + in_ljd * (!selector);
            double i2d = ijd * (selector) + iid * (!selector);
            double l2d = ljd * (selector) + lid * (!selector);
            double in_l2d = in_ljd * (selector) + in_lid * (!selector);
            
            // check if Ansatz functions on different
            // levels do not overlap and neg.
            // overlap is 1 if the functions overlap
            // if they don't overlap result is zero.
            // we mask the l2 scalar product in the end
            double q = (i1d - 1) * in_l1d;
            double p = (i1d + 1) * in_l1d;
            bool overlap = (std::max(q, (i2d - 1) * in_l2d) < std::min(p, (i2d + 1) * in_l2d));
            
            // we determine fl and fr by plugging them
            // into the sparse grids basis functions given by l2d and i2d.
            // Then we use the formular from above: 1/2 * (f_l + f_r) * 2^{-l}
            double temp_res_inner = 2.0 - fabs(l2d * q - i2d) - fabs(l2d * p - i2d);
            double temp_res_rightbound = p + q;
            double temp_res_leftbound = 2.0 - temp_res_rightbound;
            
            // now select functions evaluation depending on the
            // fact if l2d is on the boundary
            double temp_res = (temp_res_inner * (l2d != 1)) + temp_res_leftbound * ((l2d == 1) && (i2d == 0)) + temp_res_rightbound * ((l2d == 1) && (i2d == 1));
            temp_res *= (0.5 * in_l1d);
            
            double res_two = temp_res * overlap; // Now mask result
            
            // mask with lid != 1, since we have to avoid the "up" part.
            return (res_one * (lid == ljd) + res_two * (lid != ljd)) * *(lcl_q_->getPointer() + dim);
        }
        
        void OperationLTwoDotLaplaceVectorizedLinearBoundary::mult_dirichlet(sg::base::DataVector& alpha, sg::base::DataVector& result) {
            
            size_t process_i_start = all_i_start[process_index];
            size_t process_i_end = process_i_start + all_i_size[process_index];

            //std::cout << "PROCESS :" << process_index << " START :" << process_i_start << " END :" << process_i_end << " COUNT :" << (process_i_end - process_i_start) <<std::endl;

            #ifdef WITH_GLFOPS_CALCULATION
            sg::base::SGppStopwatch stopWatch;
            stopWatch.start();
            #endif


#if defined (STORE_MATRIX)
			if(! operation_result_generated_)
			{
				operation_result_generated_ = true;
				
				alpha_padded_->setAll(1.0);
#else
			std::size_t original_size = alpha.getSize();
            memcpy(alpha_padded_->getPointer(), alpha.getPointer(), original_size * sizeof(double));
#endif


#if defined(__SSE4_2__) && defined(__AVX__)

            #pragma omp parallel
            {
				std::size_t padded_size = this->level_->getNcols();               
				double* constants = this->constants_->getPointer();//{0, 0.5, 2.0 / 3.0, 1, 2};
#if ! defined (STORE_MATRIX)
				double* result_ptr_ = result_boundary_filtered_->getPointer();
#endif
				double* level_ptr_ = this->level_->getPointer();
				double* level_int_ptr_ = this->level_int_->getPointer();
				double* index_ptr_ = this->index_->getPointer();
				
				double* level_boundary_filtered_ptr_ = this->level_boundary_filtered_->getPointer();
				double* level_int_boundary_filtered_ptr_ = this->level_int_boundary_filtered_->getPointer();
				double* index_boundary_filtered_ptr_ = this->index_boundary_filtered_->getPointer();

				double* gradient_temp_ptr = gradient_temp[omp_get_thread_num()]->getPointer();
				double* l2dot_temp_ptr = l2dot_temp[omp_get_thread_num()]->getPointer();
				double* lambda_ptr_ = this->lambda_->getPointer();
				size_t temp_cols = this->storage->dim() * VECTOR_SIZE;
				
				double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
				double* lcl_q_inv_temp_ptr_ = lcl_q_inv_->getPointer();
				double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();
				
				size_t max_dims = this->storage->dim();
				
				__m256d mm_half = _mm256_broadcast_sd(constants + 1);
				__m256d mm_two_thirds = _mm256_broadcast_sd(constants + 2);
				__m256d mm_one = _mm256_broadcast_sd(constants + 3);
				__m256d mm_two = _mm256_broadcast_sd(constants + 4);
				__m256d mm_abs = _mm256_broadcast_sd(constants + 5); 

				__m256d mm_timestep_coeff = _mm256_broadcast_sd(&TimestepCoeff);
				

				size_t thr_start;
				size_t thr_end;
				sg::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end, &thr_start, &thr_end);

                for(size_t ii = thr_start; ii < thr_end; ii++)
                {
                        __m256d mm_result = _mm256_setzero_pd();
                        __m256d mm_result2 = _mm256_setzero_pd();
                        
                        double* temp_level_ptr = level_ptr_;
                        double* temp_level_int_ptr = level_int_ptr_;
                        double* temp_index_ptr = index_ptr_;
                        
                        for(size_t j = 0; j < padded_size; j+= VECTOR_SIZE * REG_BCOUNT)
                        {
#if defined (STORE_MATRIX)
						mm_result = _mm256_setzero_pd();
						mm_result2 = _mm256_setzero_pd();
#endif	
                            double* gradient_temp_ptr1 = gradient_temp_ptr;
                            double* l2dot_temp_ptr1 = l2dot_temp_ptr;
                            
                            size_t i_idx = ii * max_dims;
                            for(size_t dim = 0; dim < max_dims; dim++)
                            {
                                __m256d mm_lcl_q = _mm256_broadcast_sd(lcl_q_temp_ptr_ + dim);
                                __m256d mm_lcl_q_inv = _mm256_broadcast_sd(lcl_q_inv_temp_ptr_ + dim);
                                
								//std::cout << ii << " " << level_boundary_filtered_->getNrows() << dim << " " << level_boundary_filtered_->getNcols() << std::endl;
                                __m256d mm_lid = _mm256_broadcast_sd(level_boundary_filtered_ptr_ + i_idx);
                                __m256d mm_iid = _mm256_broadcast_sd(index_boundary_filtered_ptr_ + i_idx);
                                __m256d mm_ljd = _mm256_load_pd(temp_level_ptr);
                                __m256d mm_ijd = _mm256_load_pd(temp_index_ptr);
                                
                                __m256d mm_doGrad = _mm256_and_pd(_mm256_cmp_pd(mm_lid, mm_ljd, _CMP_EQ_OQ),
                                                                  _mm256_and_pd(_mm256_cmp_pd(mm_iid, mm_ijd, _CMP_EQ_OQ),
                                                                                _mm256_cmp_pd(mm_lid, mm_one, _CMP_NEQ_OQ))); //+5
                                
                                
                                __m256d mm_grad = _mm256_mul_pd(mm_lid, _mm256_and_pd(mm_two, mm_doGrad)); //1+1
                                
                                mm_grad = _mm256_mul_pd(mm_grad, mm_lcl_q_inv); //1
                                
                                _mm256_store_pd(gradient_temp_ptr1, mm_grad);
                                
                                
                                __m256d mm_in_lid = _mm256_broadcast_sd(level_int_boundary_filtered_ptr_ + i_idx);
                                __m256d mm_in_ljd = _mm256_load_pd(temp_level_int_ptr);
                                
                                __m256d mm_res_one = _mm256_mul_pd(mm_two_thirds, _mm256_and_pd(mm_in_lid, _mm256_and_pd(_mm256_cmp_pd(mm_iid, mm_ijd, _CMP_EQ_OQ), _mm256_cmp_pd(mm_ljd, mm_one, _CMP_NEQ_OQ)))); //1+2
                                
                                __m256d mm_selector = _mm256_cmp_pd(mm_lid, mm_ljd, _CMP_LE_OQ);//+6
                                __m256d mm_i1d = _mm256_blendv_pd(mm_iid, mm_ijd, mm_selector);
                                __m256d mm_in_l1d = _mm256_blendv_pd(mm_in_lid, mm_in_ljd, mm_selector);
                                __m256d mm_in_l2d = _mm256_blendv_pd(mm_in_ljd, mm_in_lid, mm_selector);
                                __m256d mm_i2d = _mm256_blendv_pd(mm_ijd, mm_iid, mm_selector);
                                __m256d mm_l2d = _mm256_blendv_pd(mm_ljd, mm_lid, mm_selector);
                                
                                __m256d mm_q = _mm256_mul_pd(_mm256_sub_pd(mm_i1d, mm_one), mm_in_l1d); //2 flop
                                __m256d mm_p = _mm256_mul_pd(_mm256_add_pd(mm_i1d, mm_one), mm_in_l1d); //2 flop
                                __m256d mm_overlap = _mm256_cmp_pd(_mm256_max_pd(mm_q, _mm256_mul_pd(_mm256_sub_pd(mm_i2d, mm_one), mm_in_l2d)),
                                                                   _mm256_min_pd(mm_p, _mm256_mul_pd(_mm256_add_pd(mm_i2d, mm_one), mm_in_l2d)),
                                                                   _CMP_LT_OQ); //6+1
                                
                                __m256d mm_temp_res_inner = _mm256_sub_pd(_mm256_sub_pd(mm_two,
                                                                                        _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d, mm_q), mm_i2d)))),
                                                                          _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d, mm_p), mm_i2d)))); // 8 flops
                                                                          
                                                                          __m256d mm_temp_res_rightbound = _mm256_add_pd(mm_p, mm_q); //1
                                                                          __m256d mm_temp_res_leftbound = _mm256_sub_pd(mm_two, mm_temp_res_rightbound); //1
                                                                          
                                                                          __m256d mm_temp_res = _mm256_blendv_pd(mm_temp_res_inner,
                                                                                                                 _mm256_blendv_pd(mm_temp_res_leftbound, mm_temp_res_rightbound, _mm256_cmp_pd(mm_i2d, mm_one, _CMP_EQ_OQ)),
                                                                                                                 _mm256_cmp_pd(mm_l2d, mm_one, _CMP_EQ_OQ));
                                                                          
                                                                          
                                                                          
                                                                          mm_temp_res = _mm256_mul_pd(mm_temp_res, _mm256_mul_pd(mm_half, mm_in_l1d)); // 2 flops
                                __m256d mm_res_two = _mm256_and_pd(mm_temp_res, mm_overlap); // Now mask result //+1
                                mm_selector = _mm256_cmp_pd(mm_lid, mm_ljd, _CMP_NEQ_OQ); // +1
                                
                                __m256d mm_val = _mm256_blendv_pd(mm_res_one, mm_res_two, mm_selector);  // +1
                                mm_val = _mm256_mul_pd(mm_val, mm_lcl_q); //1 flop
                                _mm256_store_pd(l2dot_temp_ptr1, mm_val);
                                
                                ////////////////////////////////////////////////////////
                                __m256d mm_ljd2 = _mm256_load_pd(temp_level_ptr + VECTOR_SIZE);
                                __m256d mm_ijd2 = _mm256_load_pd(temp_index_ptr + VECTOR_SIZE);
                                __m256d mm_doGrad2 = _mm256_and_pd(_mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_EQ_OQ),
                                                                   _mm256_and_pd(_mm256_cmp_pd(mm_iid, mm_ijd2, _CMP_EQ_OQ),
                                                                                 _mm256_cmp_pd(mm_lid, mm_one, _CMP_NEQ_OQ))); //1 // +2
                                
                                
                                __m256d mm_grad2 = _mm256_mul_pd(mm_lid, _mm256_and_pd(mm_two, mm_doGrad2)); //2
                                
                                mm_grad2 = _mm256_mul_pd(mm_grad2, mm_lcl_q_inv); //1
                                _mm256_store_pd(gradient_temp_ptr1 + temp_cols, mm_grad2);
                                
                                __m256d mm_in_ljd2 = _mm256_load_pd(temp_level_int_ptr + VECTOR_SIZE);
                                
                                __m256d mm_res_one2 = _mm256_mul_pd(mm_two_thirds, _mm256_and_pd(mm_in_lid, _mm256_and_pd(_mm256_cmp_pd(mm_iid, mm_ijd2, _CMP_EQ_OQ), _mm256_cmp_pd(mm_ljd2, mm_one, _CMP_NEQ_OQ)))); //2 // +1
                                
                                __m256d mm_selector2 = _mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_LE_OQ);
                                __m256d mm_i1d2 = _mm256_blendv_pd(mm_iid, mm_ijd2, mm_selector2);
                                __m256d mm_in_l1d2 = _mm256_blendv_pd(mm_in_lid, mm_in_ljd2, mm_selector2);
                                __m256d mm_i2d2 = _mm256_blendv_pd(mm_ijd2, mm_iid, mm_selector2);
                                __m256d mm_l2d2 = _mm256_blendv_pd(mm_ljd2, mm_lid, mm_selector2);
                                __m256d mm_in_l2d2 = _mm256_blendv_pd(mm_in_ljd2, mm_in_lid, mm_selector2);
                                
                                
                                __m256d mm_q2 = _mm256_mul_pd(_mm256_sub_pd(mm_i1d2, mm_one), mm_in_l1d2); //2 flop
                                __m256d mm_p2 = _mm256_mul_pd(_mm256_add_pd(mm_i1d2, mm_one), mm_in_l1d2); //2 flop
                                
                                __m256d mm_overlap2 = _mm256_cmp_pd(_mm256_max_pd(mm_q2, _mm256_mul_pd(_mm256_sub_pd(mm_i2d2, mm_one), mm_in_l2d2)),
                                                                    _mm256_min_pd(mm_p2, _mm256_mul_pd(_mm256_add_pd(mm_i2d2, mm_one), mm_in_l2d2)),
                                                                    _CMP_LT_OQ); //6 flop // +1
                                
                                __m256d mm_temp_res_inner2 = _mm256_sub_pd(_mm256_sub_pd(mm_two,
                                                                                         _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
                                                                           _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d2, mm_p2), mm_i2d2)))); // 8 flops
                                                                           
                                                                           __m256d mm_temp_res_rightbound2 = _mm256_add_pd(mm_p2, mm_q2);
                                                                           __m256d mm_temp_res_leftbound2 = _mm256_sub_pd(mm_two, mm_temp_res_rightbound2);
                                                                           
                                                                           __m256d mm_temp_res2 = _mm256_blendv_pd(mm_temp_res_inner2,
                                                                                                                   _mm256_blendv_pd(mm_temp_res_leftbound2, mm_temp_res_rightbound2, _mm256_cmp_pd(mm_i2d2, mm_one, _CMP_EQ_OQ)),
                                                                                                                   _mm256_cmp_pd(mm_l2d2, mm_one, _CMP_EQ_OQ));
                                                                           
                                                                           mm_temp_res2 = _mm256_mul_pd(mm_temp_res2, _mm256_mul_pd(mm_half, mm_in_l1d2)); // 2 flops
                                __m256d mm_res_two2 = _mm256_and_pd(mm_temp_res2, mm_overlap2); // Now mask result //1 flop
                                
                                mm_selector2 = _mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_NEQ_OQ);
                                
                                __m256d mm_val2 = _mm256_blendv_pd(mm_res_one2, mm_res_two2, mm_selector2);
                                
                                
                                mm_val2 = _mm256_mul_pd(mm_val2, mm_lcl_q); //1 flop
                                
                                _mm256_store_pd(l2dot_temp_ptr1 + temp_cols, mm_val2);
                                
                                
                                gradient_temp_ptr1 += VECTOR_SIZE;
                                l2dot_temp_ptr1 += VECTOR_SIZE;
                                
                                temp_level_ptr += BLOCK_LENGTH;
                                temp_level_int_ptr += BLOCK_LENGTH;
                                temp_index_ptr += BLOCK_LENGTH;
                                
                                ++i_idx;
                            }
                            
                            gradient_temp_ptr1 = gradient_temp_ptr;
							
							__m256d mm_element = _mm256_setzero_pd();
							__m256d mm_element2 = _mm256_setzero_pd();

                            
                            for(size_t d_outer = 0; d_outer < max_dims; d_outer++) // D + 2
                            {
                                mm_element = _mm256_load_pd(alpha_padded_temp_ptr_ + j);

                                mm_element2 = _mm256_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);

                                l2dot_temp_ptr1 = l2dot_temp_ptr;
                                
                                for(size_t d_inner = 0; d_inner < d_outer; d_inner++)
                                {
                                    //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
                                    __m256d mm_temp = _mm256_load_pd(l2dot_temp_ptr1);
                                    mm_element = _mm256_mul_pd(mm_element, mm_temp);
                                    __m256d mm_temp2 = _mm256_load_pd(l2dot_temp_ptr1 + temp_cols);
                                    mm_element2 = _mm256_mul_pd(mm_element2, mm_temp2);
                                    
                                    l2dot_temp_ptr1 += VECTOR_SIZE;
                                }
                                
                                l2dot_temp_ptr1 += VECTOR_SIZE;
                                
                                for(size_t d_inner = d_outer + 1; d_inner < max_dims; d_inner++)
                                {
                                    //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
                                    __m256d mm_temp = _mm256_load_pd(l2dot_temp_ptr1);
                                    mm_element = _mm256_mul_pd(mm_element, mm_temp);
                                    __m256d mm_temp2 = _mm256_load_pd(l2dot_temp_ptr1 + temp_cols);
                                    mm_element2 = _mm256_mul_pd(mm_element2, mm_temp2);
                                    
                                    l2dot_temp_ptr1 += VECTOR_SIZE;
                                }
                                
                                __m256d mm_lambda = _mm256_broadcast_sd(lambda_ptr_ + d_outer);
								
								__m256d mm_temp = _mm256_load_pd(gradient_temp_ptr1);
								__m256d mm_temp2 = _mm256_load_pd(gradient_temp_ptr1 + temp_cols);
                                
								mm_temp = _mm256_mul_pd(mm_temp, mm_lambda);
								mm_temp2 = _mm256_mul_pd(mm_temp2, mm_lambda);
								
								mm_temp = _mm256_mul_pd(mm_temp, mm_element);
								mm_temp2 = _mm256_mul_pd(mm_temp2, mm_element2);
								
								mm_temp = _mm256_mul_pd(mm_temp, mm_timestep_coeff);
								mm_temp2 = _mm256_mul_pd(mm_temp2, mm_timestep_coeff);
								
								mm_result = _mm256_add_pd(mm_result, mm_temp);
								mm_result2 = _mm256_add_pd(mm_result2, mm_temp2);
                                
                                gradient_temp_ptr1 += VECTOR_SIZE;
                            }
							
							__m256d mm_temp = _mm256_load_pd(l2dot_temp_ptr + VECTOR_SIZE * (max_dims - 1));
							__m256d mm_temp2 = _mm256_load_pd(l2dot_temp_ptr + temp_cols + VECTOR_SIZE * (max_dims - 1));
							
							mm_temp = _mm256_mul_pd(mm_element, mm_temp);
							mm_temp2 = _mm256_mul_pd(mm_element2, mm_temp2);
							
							
							mm_result = _mm256_add_pd(mm_result, mm_temp);
							mm_result2 = _mm256_add_pd(mm_result2, mm_temp2);
							
														
#if defined (STORE_MATRIX)
						double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (ii - process_i_start) * operation_result_matrix_->getNcols();
						_mm256_store_pd(operation_result_dest_ptr + j, mm_result);
						_mm256_store_pd(operation_result_dest_ptr + j + VECTOR_SIZE, mm_result2);
#endif
                        }

#if ! defined (STORE_MATRIX)
			
						mm_result = _mm256_add_pd(mm_result, mm_result2);
						mm_result = _mm256_hadd_pd(mm_result, _mm256_setzero_pd());
						
						double s_result;
						__m256d hsum = _mm256_add_pd(mm_result, _mm256_permute2f128_pd(mm_result, mm_result, 0x1));
						_mm_store_sd(&s_result, _mm_hadd_pd( _mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum) ) );
				  
						result_ptr_[ii] += s_result;
#endif
                    }
                
            }
            
#elif defined(__SSE4_2__) && !defined(__AVX__)

            #pragma omp parallel
            {
				std::size_t padded_size = this->level_->getNcols();               
				double* constants = this->constants_->getPointer();//{0, 0.5, 2.0 / 3.0, 1, 2};
#if ! defined (STORE_MATRIX)
				double* result_ptr_ = result_boundary_filtered_->getPointer();
#endif				
				double* level_ptr_ = this->level_->getPointer();
				double* level_int_ptr_ = this->level_int_->getPointer();
				double* index_ptr_ = this->index_->getPointer();
				
				double* level_boundary_filtered_ptr_ = this->level_boundary_filtered_->getPointer();
				double* level_int_boundary_filtered_ptr_ = this->level_int_boundary_filtered_->getPointer();
				double* index_boundary_filtered_ptr_ = this->index_boundary_filtered_->getPointer();

				double* gradient_temp_ptr = gradient_temp[omp_get_thread_num()]->getPointer();
				double* l2dot_temp_ptr = l2dot_temp[omp_get_thread_num()]->getPointer();
				double* lambda_ptr_ = this->lambda_->getPointer();
				size_t temp_cols = this->storage->dim() * VECTOR_SIZE;
				
				double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
				double* lcl_q_inv_temp_ptr_ = lcl_q_inv_->getPointer();
				double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();
				
				size_t max_dims = this->storage->dim();
				
				__m128d mm_half = _mm_loaddup_pd(constants + 1);
				__m128d mm_two_thirds = _mm_loaddup_pd(constants + 2);
				__m128d mm_one = _mm_loaddup_pd(constants + 3);
				__m128d mm_two = _mm_loaddup_pd(constants + 4);
				__m128d mm_abs = _mm_loaddup_pd(constants + 5); 

				__m128d mm_timestep_coeff = _mm_loaddup_pd(&TimestepCoeff);
				

				size_t thr_start;
				size_t thr_end;
				sg::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end, &thr_start, &thr_end);

                for(size_t ii = thr_start; ii < thr_end; ii++)
                {
                        __m128d mm_result = _mm_setzero_pd();
                        __m128d mm_result2 = _mm_setzero_pd();
                        
                        double* temp_level_ptr = level_ptr_;
                        double* temp_level_int_ptr = level_int_ptr_;
                        double* temp_index_ptr = index_ptr_;
                        
                        for(size_t j = 0; j < padded_size; j+= VECTOR_SIZE * REG_BCOUNT)
                        {
#if ! defined (STORE_MATRIX)
						mm_result = _mm_setzero_pd();
						mm_result2 = _mm_setzero_pd();
#endif
                            double* gradient_temp_ptr1 = gradient_temp_ptr;
                            double* l2dot_temp_ptr1 = l2dot_temp_ptr;
                            
                            size_t i_idx = ii * max_dims;
                            for(size_t dim = 0; dim < max_dims; dim++)
                            {
                                __m128d mm_lcl_q = _mm_loaddup_pd(lcl_q_temp_ptr_ + dim);
                                __m128d mm_lcl_q_inv = _mm_loaddup_pd(lcl_q_inv_temp_ptr_ + dim);
                                
								//std::cout << ii << " " << level_boundary_filtered_->getNrows() << dim << " " << level_boundary_filtered_->getNcols() << std::endl;
                                __m128d mm_lid = _mm_loaddup_pd(level_boundary_filtered_ptr_ + i_idx);
                                __m128d mm_iid = _mm_loaddup_pd(index_boundary_filtered_ptr_ + i_idx);
                                __m128d mm_ljd = _mm_load_pd(temp_level_ptr);
                                __m128d mm_ijd = _mm_load_pd(temp_index_ptr);
                                
                                __m128d mm_doGrad = _mm_and_pd(_mm_cmpeq_pd(mm_lid, mm_ljd),
                                                                  _mm_and_pd(_mm_cmpeq_pd(mm_iid, mm_ijd),
                                                                                _mm_cmpneq_pd(mm_lid, mm_one))); //+5
                                
                                
                                __m128d mm_grad = _mm_mul_pd(mm_lid, _mm_and_pd(mm_two, mm_doGrad)); //1+1
                                
                                mm_grad = _mm_mul_pd(mm_grad, mm_lcl_q_inv); //1
                                
                                _mm_store_pd(gradient_temp_ptr1, mm_grad);
                                
                                
                                __m128d mm_in_lid = _mm_loaddup_pd(level_int_boundary_filtered_ptr_ + i_idx);
                                __m128d mm_in_ljd = _mm_load_pd(temp_level_int_ptr);
                                
                                __m128d mm_res_one = _mm_mul_pd(mm_two_thirds, _mm_and_pd(mm_in_lid, _mm_and_pd(_mm_cmpeq_pd(mm_iid, mm_ijd), _mm_cmpneq_pd(mm_ljd, mm_one)))); //1+2
                                
                                __m128d mm_selector = _mm_cmple_pd(mm_lid, mm_ljd);//+6
                                __m128d mm_i1d = _mm_blendv_pd(mm_iid, mm_ijd, mm_selector);
                                __m128d mm_in_l1d = _mm_blendv_pd(mm_in_lid, mm_in_ljd, mm_selector);
                                __m128d mm_in_l2d = _mm_blendv_pd(mm_in_ljd, mm_in_lid, mm_selector);
                                __m128d mm_i2d = _mm_blendv_pd(mm_ijd, mm_iid, mm_selector);
                                __m128d mm_l2d = _mm_blendv_pd(mm_ljd, mm_lid, mm_selector);
                                
                                __m128d mm_q = _mm_mul_pd(_mm_sub_pd(mm_i1d, mm_one), mm_in_l1d); //2 flop
                                __m128d mm_p = _mm_mul_pd(_mm_add_pd(mm_i1d, mm_one), mm_in_l1d); //2 flop
                                __m128d mm_overlap = _mm_cmplt_pd(_mm_max_pd(mm_q, _mm_mul_pd(_mm_sub_pd(mm_i2d, mm_one), mm_in_l2d)),
                                                                   _mm_min_pd(mm_p, _mm_mul_pd(_mm_add_pd(mm_i2d, mm_one), mm_in_l2d))); //6+1
                                
                                __m128d mm_temp_res_inner = _mm_sub_pd(_mm_sub_pd(mm_two,
                                                                                        _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d, mm_q), mm_i2d)))),
                                                                          _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d, mm_p), mm_i2d)))); // 8 flops
                                                                          
                                                                          __m128d mm_temp_res_rightbound = _mm_add_pd(mm_p, mm_q); //1
                                                                          __m128d mm_temp_res_leftbound = _mm_sub_pd(mm_two, mm_temp_res_rightbound); //1
                                                                          
                                                                          __m128d mm_temp_res = _mm_blendv_pd(mm_temp_res_inner,
                                                                                                                 _mm_blendv_pd(mm_temp_res_leftbound, mm_temp_res_rightbound, _mm_cmpeq_pd(mm_i2d, mm_one)),
                                                                                                                 _mm_cmpeq_pd(mm_l2d, mm_one));
                                                                          
                                                                          
                                                                          
                                                                          mm_temp_res = _mm_mul_pd(mm_temp_res, _mm_mul_pd(mm_half, mm_in_l1d)); // 2 flops
                                __m128d mm_res_two = _mm_and_pd(mm_temp_res, mm_overlap); // Now mask result //+1
                                mm_selector = _mm_cmpneq_pd(mm_lid, mm_ljd); // +1
                                
                                __m128d mm_val = _mm_blendv_pd(mm_res_one, mm_res_two, mm_selector);  // +1
                                mm_val = _mm_mul_pd(mm_val, mm_lcl_q); //1 flop
                                _mm_store_pd(l2dot_temp_ptr1, mm_val);
                                
                                ////////////////////////////////////////////////////////
                                __m128d mm_ljd2 = _mm_load_pd(temp_level_ptr + VECTOR_SIZE);
                                __m128d mm_ijd2 = _mm_load_pd(temp_index_ptr + VECTOR_SIZE);
                                __m128d mm_doGrad2 = _mm_and_pd(_mm_cmpeq_pd(mm_lid, mm_ljd2),
                                                                   _mm_and_pd(_mm_cmpeq_pd(mm_iid, mm_ijd2),
                                                                                 _mm_cmpneq_pd(mm_lid, mm_one))); //1 // +2
                                
                                
                                __m128d mm_grad2 = _mm_mul_pd(mm_lid, _mm_and_pd(mm_two, mm_doGrad2)); //2
                                
                                mm_grad2 = _mm_mul_pd(mm_grad2, mm_lcl_q_inv); //1
                                _mm_store_pd(gradient_temp_ptr1 + temp_cols, mm_grad2);
                                
                                __m128d mm_in_ljd2 = _mm_load_pd(temp_level_int_ptr + VECTOR_SIZE);
                                
                                __m128d mm_res_one2 = _mm_mul_pd(mm_two_thirds, _mm_and_pd(mm_in_lid, _mm_and_pd(_mm_cmpeq_pd(mm_iid, mm_ijd2), _mm_cmpneq_pd(mm_ljd2, mm_one)))); //2 // +1
                                
                                __m128d mm_selector2 = _mm_cmple_pd(mm_lid, mm_ljd2);
                                __m128d mm_i1d2 = _mm_blendv_pd(mm_iid, mm_ijd2, mm_selector2);
                                __m128d mm_in_l1d2 = _mm_blendv_pd(mm_in_lid, mm_in_ljd2, mm_selector2);
                                __m128d mm_i2d2 = _mm_blendv_pd(mm_ijd2, mm_iid, mm_selector2);
                                __m128d mm_l2d2 = _mm_blendv_pd(mm_ljd2, mm_lid, mm_selector2);
                                __m128d mm_in_l2d2 = _mm_blendv_pd(mm_in_ljd2, mm_in_lid, mm_selector2);
                                
                                
                                __m128d mm_q2 = _mm_mul_pd(_mm_sub_pd(mm_i1d2, mm_one), mm_in_l1d2); //2 flop
                                __m128d mm_p2 = _mm_mul_pd(_mm_add_pd(mm_i1d2, mm_one), mm_in_l1d2); //2 flop
                                
                                __m128d mm_overlap2 = _mm_cmplt_pd(_mm_max_pd(mm_q2, _mm_mul_pd(_mm_sub_pd(mm_i2d2, mm_one), mm_in_l2d2)),
                                                                    _mm_min_pd(mm_p2, _mm_mul_pd(_mm_add_pd(mm_i2d2, mm_one), mm_in_l2d2))); //6 flop // +1
                                
                                __m128d mm_temp_res_inner2 = _mm_sub_pd(_mm_sub_pd(mm_two,
                                                                                         _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
                                                                           _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d2, mm_p2), mm_i2d2)))); // 8 flops
                                                                           
                                                                           __m128d mm_temp_res_rightbound2 = _mm_add_pd(mm_p2, mm_q2);
                                                                           __m128d mm_temp_res_leftbound2 = _mm_sub_pd(mm_two, mm_temp_res_rightbound2);
                                                                           
                                                                           __m128d mm_temp_res2 = _mm_blendv_pd(mm_temp_res_inner2,
                                                                                                                   _mm_blendv_pd(mm_temp_res_leftbound2, mm_temp_res_rightbound2, _mm_cmpeq_pd(mm_i2d2, mm_one)),
                                                                                                                   _mm_cmpeq_pd(mm_l2d2, mm_one));
                                                                           
                                                                           mm_temp_res2 = _mm_mul_pd(mm_temp_res2, _mm_mul_pd(mm_half, mm_in_l1d2)); // 2 flops
                                __m128d mm_res_two2 = _mm_and_pd(mm_temp_res2, mm_overlap2); // Now mask result //1 flop
                                
                                mm_selector2 = _mm_cmpneq_pd(mm_lid, mm_ljd2);
                                
                                __m128d mm_val2 = _mm_blendv_pd(mm_res_one2, mm_res_two2, mm_selector2);
                                
                                
                                mm_val2 = _mm_mul_pd(mm_val2, mm_lcl_q); //1 flop
                                
                                _mm_store_pd(l2dot_temp_ptr1 + temp_cols, mm_val2);
                                
                                
                                gradient_temp_ptr1 += VECTOR_SIZE;
                                l2dot_temp_ptr1 += VECTOR_SIZE;
                                
                                temp_level_ptr += BLOCK_LENGTH;
                                temp_level_int_ptr += BLOCK_LENGTH;
                                temp_index_ptr += BLOCK_LENGTH;
                                
                                ++i_idx;
                            }
                            
                            gradient_temp_ptr1 = gradient_temp_ptr;
							
							__m128d mm_element = _mm_setzero_pd();
							__m128d mm_element2 = _mm_setzero_pd();

                            
                            for(size_t d_outer = 0; d_outer < max_dims; d_outer++) // D + 2
                            {
                                mm_element = _mm_load_pd(alpha_padded_temp_ptr_ + j);

                                mm_element2 = _mm_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);

                                l2dot_temp_ptr1 = l2dot_temp_ptr;
                                
                                for(size_t d_inner = 0; d_inner < d_outer; d_inner++)
                                {
                                    //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
                                    __m128d mm_temp = _mm_load_pd(l2dot_temp_ptr1);
                                    mm_element = _mm_mul_pd(mm_element, mm_temp);
                                    __m128d mm_temp2 = _mm_load_pd(l2dot_temp_ptr1 + temp_cols);
                                    mm_element2 = _mm_mul_pd(mm_element2, mm_temp2);
                                    
                                    l2dot_temp_ptr1 += VECTOR_SIZE;
                                }
                                
                                l2dot_temp_ptr1 += VECTOR_SIZE;
                                
                                for(size_t d_inner = d_outer + 1; d_inner < max_dims; d_inner++)
                                {
                                    //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
                                    __m128d mm_temp = _mm_load_pd(l2dot_temp_ptr1);
                                    mm_element = _mm_mul_pd(mm_element, mm_temp);
                                    __m128d mm_temp2 = _mm_load_pd(l2dot_temp_ptr1 + temp_cols);
                                    mm_element2 = _mm_mul_pd(mm_element2, mm_temp2);
                                    
                                    l2dot_temp_ptr1 += VECTOR_SIZE;
                                }
                                
                                __m128d mm_lambda = _mm_loaddup_pd(lambda_ptr_ + d_outer);
								
								__m128d mm_temp = _mm_load_pd(gradient_temp_ptr1);
								__m128d mm_temp2 = _mm_load_pd(gradient_temp_ptr1 + temp_cols);
                                
								mm_temp = _mm_mul_pd(mm_temp, mm_lambda);
								mm_temp2 = _mm_mul_pd(mm_temp2, mm_lambda);
								
								mm_temp = _mm_mul_pd(mm_temp, mm_element);
								mm_temp2 = _mm_mul_pd(mm_temp2, mm_element2);
								
								mm_temp = _mm_mul_pd(mm_temp, mm_timestep_coeff);
								mm_temp2 = _mm_mul_pd(mm_temp2, mm_timestep_coeff);
								
								mm_result = _mm_add_pd(mm_result, mm_temp);
								mm_result2 = _mm_add_pd(mm_result2, mm_temp2);
                                
                                gradient_temp_ptr1 += VECTOR_SIZE;
                            }
							__m128d mm_temp = _mm_load_pd(l2dot_temp_ptr + VECTOR_SIZE * (max_dims - 1));
							__m128d mm_temp2 = _mm_load_pd(l2dot_temp_ptr + temp_cols + VECTOR_SIZE * (max_dims - 1));
							
							mm_temp = _mm_mul_pd(mm_element, mm_temp);
							mm_temp2 = _mm_mul_pd(mm_element2, mm_temp2);
							
							
							mm_result = _mm_add_pd(mm_result, mm_temp);
							mm_result2 = _mm_add_pd(mm_result2, mm_temp2);
							
#if defined (STORE_MATRIX)
							double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (ii - process_i_start) * operation_result_matrix_->getNcols();
							_mm_store_pd(operation_result_dest_ptr + j, mm_result);
							_mm_store_pd(operation_result_dest_ptr + j + VECTOR_SIZE, mm_result2);
#endif
                        }

#if ! defined (STORE_MATRIX)                        
						mm_result = _mm_add_pd(mm_result, mm_result2);
						mm_result = _mm_hadd_pd(mm_result, _mm_setzero_pd());
						
						double s_result = 0.0;
						_mm_store_sd(&s_result, mm_result);

						result_ptr_[ii] += s_result;
#endif							
				}
			}
#else
#error "Needs SSE4.2 or AVX to compile"
#endif


#if defined (STORE_MATRIX)
	}
	
	std::size_t original_size = alpha.getSize();
	memcpy(alpha_padded_->getPointer(), alpha.getPointer(), original_size * sizeof(double));
	
#pragma omp parallel
{
		size_t padded_size = this->operation_result_matrix_->getNcols();
		size_t thr_start;
		size_t thr_end;
		sg::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end, &thr_start, &thr_end);
		
		
        double* alpha_padded_ptr_ = alpha_padded_->getPointer();
		double* result_boundary_filtered_ptr = result_boundary_filtered_->getPointer();
		
		for(size_t ii = thr_start; ii < thr_end; ii++)
		{
			double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (ii - process_i_start) * operation_result_matrix_->getNcols();
			
#if defined(__SSE4_2__) && defined(__AVX__)
            __m256d mm_element = _mm256_setzero_pd();
						
			for(size_t j = 0; j < padded_size; j += VECTOR_SIZE)
			{
				__m256d mm_temp1 = _mm256_load_pd(alpha_padded_ptr_ + j);
				
				__m256d mm_temp2 = _mm256_load_pd(operation_result_dest_ptr + j);
				
				mm_element = _mm256_add_pd(mm_element, _mm256_mul_pd(mm_temp1, mm_temp2));
			}
			
			 __m256d hsum = _mm256_add_pd(mm_element, _mm256_permute2f128_pd(mm_element, mm_element, 0x1));
			  _mm_store_sd(result_boundary_filtered_ptr + ii, _mm_hadd_pd( _mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum) ) );
			
#elif defined(__SSE4_2__) && !defined(__AVX__)
            __m128d mm_element = _mm_setzero_pd();
						
			for(size_t j = 0; j < padded_size; j += VECTOR_SIZE)
			{
				__m128d mm_temp1 = _mm_load_pd(alpha_padded_ptr_ + j);
				
				__m128d mm_temp2 = _mm_load_pd(operation_result_dest_ptr + j);
				
				mm_element = _mm_add_pd(mm_element, _mm_mul_pd(mm_temp1, mm_temp2));
			}
			
			mm_element = _mm_hadd_pd(mm_element, _mm_setzero_pd());
			_mm_store_sd(result_boundary_filtered_ptr + ii, mm_element);
#else
			double element = 0.0;
			
			for(size_t j = 0 ;j < storage->size() ; ++j)
			{
				element += alpha[j] * *(operation_result_dest_ptr + j);
			}
			
			result_boundary_filtered_ptr[ii] = element;
#endif

		}
	}
#endif

	double* result_boundary_filtered_ptr = result_boundary_filtered_->getPointer();
#ifdef USE_MPI            
	/*
	 MPI_Alltoallv(result_boundary_filtered_ptr, send_size.data(), send_start.data(), MPI_DOUBLE,
				   result_boundary_filtered_ptr, recv_size.data(), recv_start.data(), MPI_DOUBLE,
				   MPI_COMM_WORLD);
	*/						   
	MPI_Allgatherv(MPI_IN_PLACE, send_size[0], MPI_DOUBLE,
				   result_boundary_filtered_ptr, recv_size.data(), recv_start.data(), 
				   MPI_DOUBLE, MPI_COMM_WORLD);
#endif

	//move to result vector
	double* result_ptr = result.getPointer();
	
	#pragma omp parallel for
	for(size_t i = 0; i < i_boundary_filtered.size(); ++i)
	{
		result_ptr[i_boundary_filtered[i]] = result_boundary_filtered_ptr[i];
	}
	
#ifdef WITH_GLFOPS_CALCULATION
	double needed_time = stopWatch.stop();
	
	size_t flop = ((28) * storage->dim() + storage->dim() * storage->dim()) * storage->size() * all_i_size[process_index];
	double gflops = (flop / needed_time) / 1000000000;
	size_t gop = ((54) * storage->dim() + storage->dim() * storage->dim()) * storage->size() * all_i_size[process_index];
	double gops = (gop / needed_time) / 1000000000;
	
	for(int i = 0; i < process_count ; i++)
	{
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if(i == process_index)
		{            
			//if (process_index == 0)
			std::cout << "[PROCESS :" << process_index << "] GFLOPS :" << gflops << " = (" << (flop/1000000000) << " / "<< needed_time << ")  GOPS :" << gops  << std::endl;
		}
		std::cout.flush();
	}
#endif

}

void OperationLTwoDotLaplaceVectorizedLinearBoundary::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
    result.setAll(0.0);
	result_boundary_filtered_->setAll(0.0);
    bool dirichlet = true;
    
    // fill q array
    for (size_t d = 0; d < this->storage->dim(); d++) {
        sg::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
        *(lcl_q_->getPointer() + d) = boundingBox->getIntervalWidth(d);
        *(lcl_q_inv_->getPointer() + d) = 1.0 / boundingBox->getIntervalWidth(d);
        dirichlet = dirichlet && boundingBox->hasDirichletBoundaryLeft(d);
        dirichlet = dirichlet && boundingBox->hasDirichletBoundaryRight(d);
    }
    
    if (dirichlet) {
        mult_dirichlet(alpha, result);
    } else {
        throw new sg::base::operation_exception("OperationLTwoDotLaplaceVectorizedLinearBoundary::mult : This method is only available on grids with Dirichlet boundaries in all dimensions!");
    }
}

}

}
