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

#include "parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotProductVectorizedLinear.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/tools/TypesParallel.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "base/exception/operation_exception.hpp"

#include <cmath>

#include <cstring>
#include <iostream>

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

        OperationLTwoDotProductVectorizedLinear::OperationLTwoDotProductVectorizedLinear(sg::base::GridStorage* storage)
			: storage(storage),
			 level_(NULL),
			 level_int_(NULL),
			 index_(NULL),
			 lcl_q_(NULL),
			 alpha_padded_(NULL),
			 constants_(NULL),
			 gradient_temp(NULL),
			 l2dot_temp(NULL)
#if defined(STORE_MATRIX)
			,
			operation_result_matrix_(NULL),
			operation_result_generated_(false)
#endif 
		{
std::cout<<"IN CONSTRUSTOR: OperationLTwoDotProductVectorizedLinear" << std::endl;
            init_constants();
			init_grid_storage();
        }

        OperationLTwoDotProductVectorizedLinear::~OperationLTwoDotProductVectorizedLinear() {
            delete this->level_;
            delete this->level_int_;
            delete this->index_;
            delete lcl_q_;

            delete this->alpha_padded_;
			delete this->constants_;
			
#pragma omp parallel
{
			delete gradient_temp[omp_get_thread_num()];
			delete l2dot_temp[omp_get_thread_num()];
}
            delete gradient_temp;
			delete l2dot_temp;
			
#if defined (STORE_MATRIX)			
			delete operation_result_matrix_;
#endif
        }
        
        void OperationLTwoDotProductVectorizedLinear::reset() {
			init_grid_storage();
		}
        
        void OperationLTwoDotProductVectorizedLinear::init_constants() {
		
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

        void OperationLTwoDotProductVectorizedLinear::init_grid_storage() {
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
			lcl_q_ = new sg::base::DataVector(this->storage->dim());
            double* lcl_q_ptr_ = lcl_q_->getPointer();
			
            
            storage->getLevelIndexArraysForEvalTLBOptimized(*(this->level_), *(this->index_), sg::parallel::X86SIMD, BLOCK_LENGTH);
            storage->getLevelForIntegralTLBOptimized(*(this->level_int_), sg::parallel::X86SIMD, BLOCK_LENGTH);
			
            std::size_t padded_size = this->level_->getNcols();
			if(this->alpha_padded_)
				delete this->alpha_padded_;
            this->alpha_padded_ = new sg::base::DataVector(padded_size);
            this->alpha_padded_->setAll(0.0);
            
            
            size_t single_process_portion = (this->storage->size() / process_count) + 1;
            
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
                int process_portion = (i == process_count - 1)? (int) (this->storage->size() - i * single_process_portion)
																: (int) (single_process_portion);
                
				process_portion = std::min<int>(process_portion, (int) (this->storage->size() - i * single_process_portion));
				process_portion = std::max<int>(process_portion, 0);
				
				all_i_size.push_back(process_portion);
            }
            
            for(int i = 0; i < process_count; ++i)
            {
                send_start.push_back((int) (single_process_portion * process_index));
                int process_send_size = (process_index == process_count -1)?  (int) (this->storage->size() - single_process_portion * process_index)
																			: (int) (single_process_portion);
																			
                process_send_size = std::min<int>(process_send_size, (int) (this->storage->size() - single_process_portion * process_index));
				process_send_size = std::max<int>(process_send_size, 0);
				//process_send_size = (i == process_index)? 0: process_send_size;
				
                send_size.push_back(process_send_size);
                
                recv_start.push_back((int) (single_process_portion * i));
                int process_recv_size = (i == process_count -1)?  (int) (this->storage->size() - single_process_portion * i)
																: (int) (single_process_portion);

                process_recv_size = std::min<int>(process_recv_size, (int) (this->storage->size() - single_process_portion * i));
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
            }
			
#if defined (STORE_MATRIX)
			size_t result_matrix_rows = all_i_size[process_index];
			size_t result_matrix_cols = this->level_->getNcols();
			
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
				sprintf(exception_string, "OperationLTwoDotProductVectorizedLinear::init : More memory (= %i GB) needed to store the operation matrix, Please set the SGPP_PDE_MATRIX_SIZE_GB environment variable!", (int) matrix_needed_size_gb);
				
				std::cerr << exception_string << std::endl;
		        throw new sg::base::operation_exception(exception_string);
			}
			
			
			if(operation_result_matrix_)
				delete operation_result_matrix_;
				
			operation_result_matrix_ = new sg::base::DataMatrix(result_matrix_rows, result_matrix_cols);
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

        double OperationLTwoDotProductVectorizedLinear::l2dot(size_t i, size_t j, size_t dim) {
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
            double res_one = (2.0 / 3.0) * in_lid * (iid == ijd);

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

            #define ALEX
            #ifdef BENNI
            // We determine the distance and mirrow the
            // it to the left, the we apply the gradient of the
            // linear overlapped part. Finally we transform
            // back the our actual basis function by
            // multiplying with 2^{-l}
            double diff = (i1d * in_l1d) - (i2d * in_l2d);
            double temp_res = fabs(diff - in_l1d) + fabs(diff + in_l1d) - fabs(diff);
            temp_res *= l2d;
            temp_res = (1 - temp_res) * in_l1d;
            #endif
            #ifdef ALEX
            // we determine fl and fr by plugging them
            // into the sparse grids basis functions given by l2d and i2d.
            // Then we use the formular from above: 1/2 * (f_l + f_r) * 2^{-l}
            double temp_res = 2.0 - fabs(l2d * q - i2d) - fabs(l2d * p - i2d);
            temp_res *= (0.5 * in_l1d);
            #endif
            double res_two = temp_res * overlap; // Now mask result

            return (res_one * (lid == ljd) + res_two * (lid != ljd)) * *(lcl_q_->getPointer() + dim);
        }

        void OperationLTwoDotProductVectorizedLinear::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
            result.setAll(0.0);

            size_t process_i_start = all_i_start[process_index];
            size_t process_i_end = process_i_start + all_i_size[process_index];

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
				double* result_ptr_ = result.getPointer();
#endif
				double* level_ptr_ = this->level_->getPointer();
				double* level_int_ptr_ = this->level_int_->getPointer();
				double* index_ptr_ = this->index_->getPointer();
				double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
				double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();

                size_t max_dims = this->storage->dim();
                size_t page_cap_rounded = max_dims * BLOCK_LENGTH;
                
				__m256d mm_half = _mm256_broadcast_sd(constants + 1);
				__m256d mm_two_thirds = _mm256_broadcast_sd(constants + 2);
				__m256d mm_one = _mm256_broadcast_sd(constants + 3);
				__m256d mm_two = _mm256_broadcast_sd(constants + 4);
				__m256d mm_abs = _mm256_broadcast_sd(constants + 5);
				


				size_t thr_start;
				size_t thr_end;
				sg::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end, &thr_start, &thr_end);

                for(size_t i = thr_start; i < thr_end; i++)
                {
                    size_t i_page = i / BLOCK_LENGTH;
                    size_t i_offset = i_page *  page_cap_rounded + i % BLOCK_LENGTH;

                    double* temp_level_ptr = level_ptr_;
                    double* temp_level_int_ptr = level_int_ptr_;
                    double* temp_index_ptr = index_ptr_;

                    for(size_t j = 0; j < padded_size; j+= VECTOR_SIZE * REG_BCOUNT)
                    {
                        size_t i_idx = i_offset;

                        __m256d mm_element = _mm256_load_pd(alpha_padded_temp_ptr_ + j);
                        __m256d mm_element2 = _mm256_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);

                        for(size_t dim = 0; dim < max_dims; dim++)
                        {
                            __m256d mm_lcl_q = _mm256_broadcast_sd(lcl_q_temp_ptr_ + dim);

                            __m256d mm_lid = _mm256_broadcast_sd(level_ptr_ + i_idx);
                            __m256d mm_iid = _mm256_broadcast_sd(index_ptr_ + i_idx);
                            __m256d mm_ljd = _mm256_load_pd(temp_level_ptr);
                            __m256d mm_ijd = _mm256_load_pd(temp_index_ptr);

                            __m256d mm_in_lid = _mm256_broadcast_sd(level_int_ptr_ + i_idx);
                            __m256d mm_in_ljd = _mm256_load_pd(temp_level_int_ptr);

                            __m256d mm_res_one = _mm256_mul_pd(mm_two_thirds, _mm256_and_pd(mm_in_lid, _mm256_cmp_pd(mm_iid, mm_ijd, _CMP_EQ_OQ))); //1+2

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

                            __m256d mm_temp_res = _mm256_sub_pd(_mm256_sub_pd(mm_two,
                                                                            _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d, mm_q), mm_i2d)))),
                                                                _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d, mm_p), mm_i2d)))); // 8 flops
                            mm_temp_res = _mm256_mul_pd(mm_temp_res, _mm256_mul_pd(mm_half, mm_in_l1d)); // 2 flops
                            __m256d mm_res_two = _mm256_and_pd(mm_temp_res, mm_overlap); // Now mask result //+1
                            mm_selector = _mm256_cmp_pd(mm_lid, mm_ljd, _CMP_NEQ_OQ); // +1

                            __m256d mm_val = _mm256_blendv_pd(mm_res_one, mm_res_two, mm_selector);  // +1
                            mm_val = _mm256_mul_pd(mm_val, mm_lcl_q); //1 flop

                            mm_element = _mm256_mul_pd(mm_element, mm_val);// 1 flop

                            ////////////////////////////////////////////////////////
                            __m256d mm_ljd2 = _mm256_load_pd(temp_level_ptr + VECTOR_SIZE);
                            __m256d mm_ijd2 = _mm256_load_pd(temp_index_ptr + VECTOR_SIZE);

                            __m256d mm_in_ljd2 = _mm256_load_pd(temp_level_int_ptr + VECTOR_SIZE);

                            __m256d mm_res_one2 = _mm256_mul_pd(mm_two_thirds, _mm256_and_pd(mm_in_lid, _mm256_cmp_pd(mm_iid, mm_ijd2, _CMP_EQ_OQ))); //2 // +1
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

                            __m256d mm_temp_res2 = _mm256_sub_pd(_mm256_sub_pd(mm_two,
                                                                                _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
                                                                _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d2, mm_p2), mm_i2d2)))); // 8 flops

                            mm_temp_res2 = _mm256_mul_pd(mm_temp_res2, _mm256_mul_pd(mm_half, mm_in_l1d2)); // 2 flops
                            __m256d mm_res_two2 = _mm256_and_pd(mm_temp_res2, mm_overlap2); // Now mask result //1 flop

                            mm_selector2 = _mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_NEQ_OQ);

                            __m256d mm_val2 = _mm256_blendv_pd(mm_res_one2, mm_res_two2, mm_selector2);
                            mm_val2 = _mm256_mul_pd(mm_val2, mm_lcl_q); //1 flop

                            mm_element2 = _mm256_mul_pd(mm_element2, mm_val2);

                            temp_level_ptr += BLOCK_LENGTH;
                            temp_level_int_ptr += BLOCK_LENGTH;
                            temp_index_ptr += BLOCK_LENGTH;

                            i_idx += BLOCK_LENGTH;
                        }
#if defined (STORE_MATRIX)
						double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (i - process_i_start) * operation_result_matrix_->getNcols();
						_mm256_store_pd(operation_result_dest_ptr + j, mm_element);
						_mm256_store_pd(operation_result_dest_ptr + j + VECTOR_SIZE, mm_element2);
#else
                        __m256d mm_result = mm_element;
                        mm_result = _mm256_add_pd(mm_result, mm_element2);

                        mm_result = _mm256_hadd_pd(mm_result, _mm256_setzero_pd());
					
						double s_result;
						__m256d hsum = _mm256_add_pd(mm_result, _mm256_permute2f128_pd(mm_result, mm_result, 0x1));
						_mm_store_sd(&s_result, _mm_hadd_pd( _mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum) ) );
				  
						result_ptr_[i] += s_result;
#endif						
                    }
                }
            }

#elif defined(__SSE4_2__) && !defined(__AVX__)

            #pragma omp parallel
            {
				std::size_t padded_size = this->level_->getNcols();               
				double* constants = this->constants_->getPointer();//{0, 0.5, 2.0 / 3.0, 1, 2};
#if ! defined (STORE_MATRIX)
				double* result_ptr_ = result.getPointer();
#endif
				double* level_ptr_ = this->level_->getPointer();
				double* level_int_ptr_ = this->level_int_->getPointer();
				double* index_ptr_ = this->index_->getPointer();

				double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
				double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();

                size_t max_dims = this->storage->dim();
                size_t page_cap_rounded = max_dims * BLOCK_LENGTH;
                
				__m128d mm_half = _mm_loaddup_pd(constants + 1);
				__m128d mm_two_thirds = _mm_loaddup_pd(constants + 2);
				__m128d mm_one = _mm_loaddup_pd(constants + 3);
				__m128d mm_two = _mm_loaddup_pd(constants + 4);
				__m128d mm_abs = _mm_loaddup_pd(constants + 5);
				


				size_t thr_start;
				size_t thr_end;
				sg::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end, &thr_start, &thr_end);

                for(size_t i = thr_start; i < thr_end; i++)
                {
                    size_t i_page = i / BLOCK_LENGTH;
                    size_t i_offset = i_page *  page_cap_rounded + i % BLOCK_LENGTH;

                    double* temp_level_ptr = level_ptr_;
                    double* temp_level_int_ptr = level_int_ptr_;
                    double* temp_index_ptr = index_ptr_;

                    for(size_t j = 0; j < padded_size; j+= VECTOR_SIZE * REG_BCOUNT)
                    {
                        size_t i_idx = i_offset;

                        __m128d mm_element = _mm_load_pd(alpha_padded_temp_ptr_ + j);
                        __m128d mm_element2 = _mm_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);

                        for(size_t dim = 0; dim < max_dims; dim++)
                        {
                            __m128d mm_lcl_q = _mm_loaddup_pd(lcl_q_temp_ptr_ + dim);

                            __m128d mm_lid = _mm_loaddup_pd(level_ptr_ + i_idx);
                            __m128d mm_iid = _mm_loaddup_pd(index_ptr_ + i_idx);
                            __m128d mm_ljd = _mm_load_pd(temp_level_ptr);
                            __m128d mm_ijd = _mm_load_pd(temp_index_ptr);

                            __m128d mm_in_lid = _mm_loaddup_pd(level_int_ptr_ + i_idx);
                            __m128d mm_in_ljd = _mm_load_pd(temp_level_int_ptr);

                            __m128d mm_res_one = _mm_mul_pd(mm_two_thirds, _mm_and_pd(mm_in_lid, _mm_cmpeq_pd(mm_iid, mm_ijd))); //1+2

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

                            __m128d mm_temp_res = _mm_sub_pd(_mm_sub_pd(mm_two,
                                                                        _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d, mm_q), mm_i2d)))),
                                                             _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d, mm_p), mm_i2d)))); // 8 flops
                            mm_temp_res = _mm_mul_pd(mm_temp_res, _mm_mul_pd(mm_half, mm_in_l1d)); // 2 flops
                            __m128d mm_res_two = _mm_and_pd(mm_temp_res, mm_overlap); // Now mask result //+1
                            mm_selector = _mm_cmpneq_pd(mm_lid, mm_ljd); // +1

                            __m128d mm_val = _mm_blendv_pd(mm_res_one, mm_res_two, mm_selector);  // +1
                            mm_val = _mm_mul_pd(mm_val, mm_lcl_q); //1 flop
                            mm_element = _mm_mul_pd(mm_element, mm_val);

                            ////////////////////////////////////////////////////////
                            __m128d mm_ljd2 = _mm_load_pd(temp_level_ptr + VECTOR_SIZE);
                            __m128d mm_ijd2 = _mm_load_pd(temp_index_ptr + VECTOR_SIZE);


                            __m128d mm_in_ljd2 = _mm_load_pd(temp_level_int_ptr + VECTOR_SIZE);

                            __m128d mm_res_one2 = _mm_mul_pd(mm_two_thirds, _mm_and_pd(mm_in_lid, _mm_cmpeq_pd(mm_iid, mm_ijd2))); //2 // +1
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

                            __m128d mm_temp_res2 = _mm_sub_pd(_mm_sub_pd(mm_two,
                                                                        _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
                                                                _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d2, mm_p2), mm_i2d2)))); // 8 flops


                            mm_temp_res2 = _mm_mul_pd(mm_temp_res2, _mm_mul_pd(mm_half, mm_in_l1d2)); // 2 flops
                            __m128d mm_res_two2 = _mm_and_pd(mm_temp_res2, mm_overlap2); // Now mask result //1 flop

                            mm_selector2 = _mm_cmpneq_pd(mm_lid, mm_ljd2);

                            __m128d mm_val2 = _mm_blendv_pd(mm_res_one2, mm_res_two2, mm_selector2);
                            mm_val2 = _mm_mul_pd(mm_val2, mm_lcl_q); //1 flop
                            mm_element2 = _mm_mul_pd(mm_element2, mm_val2);

                            temp_level_ptr += BLOCK_LENGTH;
                            temp_level_int_ptr += BLOCK_LENGTH;
                            temp_index_ptr += BLOCK_LENGTH;

                            i_idx += BLOCK_LENGTH;
                        }
#if defined (STORE_MATRIX)
						double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (i - process_i_start) * operation_result_matrix_->getNcols();
						_mm_store_pd(operation_result_dest_ptr + j, mm_element);
						_mm_store_pd(operation_result_dest_ptr + j + VECTOR_SIZE, mm_element2);
#else
                         __m128d mm_result = mm_element;
                        mm_result = _mm_add_pd(mm_result, mm_element2);

                        mm_result = _mm_hadd_pd(mm_result, _mm_setzero_pd());

						double s_result = 0.0;
						_mm_store_sd(&s_result, mm_result);

						result_ptr_[i] += s_result;
#endif
                    }
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
		double* result_ptr = result.getPointer();
		size_t padded_size = this->operation_result_matrix_->getNcols();
		size_t thr_start;
		size_t thr_end;
		sg::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end, &thr_start, &thr_end);
		
		
        double* alpha_padded_ptr_ = alpha_padded_->getPointer();
		
		for(size_t i = thr_start; i < thr_end; i++)
		{
			double* operation_result_dest_ptr = operation_result_matrix_->getPointer() + (i - process_i_start) * operation_result_matrix_->getNcols();
			
#if defined(__SSE4_2__) && defined(__AVX__)
            __m256d mm_element = _mm256_setzero_pd();
						
			for(size_t j = 0; j < padded_size; j += VECTOR_SIZE)
			{
				__m256d mm_temp1 = _mm256_load_pd(alpha_padded_ptr_ + j);
				
				__m256d mm_temp2 = _mm256_load_pd(operation_result_dest_ptr + j);
				
				mm_element = _mm256_add_pd(mm_element, _mm256_mul_pd(mm_temp1, mm_temp2));
			}
			
			 __m256d hsum = _mm256_add_pd(mm_element, _mm256_permute2f128_pd(mm_element, mm_element, 0x1));
			  _mm_store_sd(result_ptr + i, _mm_hadd_pd( _mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum) ) );
			
#elif defined(__SSE4_2__) && !defined(__AVX__)
            __m128d mm_element = _mm_setzero_pd();
						
			for(size_t j = 0; j < padded_size; j += VECTOR_SIZE)
			{
				__m128d mm_temp1 = _mm_load_pd(alpha_padded_ptr_ + j);
				
				__m128d mm_temp2 = _mm_load_pd(operation_result_dest_ptr + j);
				
				mm_element = _mm_add_pd(mm_element, _mm_mul_pd(mm_temp1, mm_temp2));
			}
			
			mm_element = _mm_hadd_pd(mm_element, _mm_setzero_pd());
			_mm_store_sd(result_ptr + i, mm_element);
#else
			double element = 0.0;
			
			for(size_t j = 0 ;j < storage->size() ; ++j)
			{
				element += alpha[j] * *(operation_result_dest_ptr + j);
			}
			
			result[i] = element;
#endif

		}
	}
#endif
#ifdef USE_MPI
			double* result_ptr = result.getPointer();
	/*
	MPI_Alltoallv(result_ptr, send_size.data(), send_start.data(), MPI_DOUBLE,
				   result_ptr, recv_size.data(), recv_start.data(), MPI_DOUBLE,
				   MPI_COMM_WORLD);
	*/
	MPI_Allgatherv(MPI_IN_PLACE, send_size[0], MPI_DOUBLE,
				   result_ptr, recv_size.data(), recv_start.data(), 
				   MPI_DOUBLE, MPI_COMM_WORLD);
#endif


#ifdef WITH_GLFOPS_CALCULATION
	size_t flop = (23 * storage->dim() ) * storage->size() * storage->size();
	size_t gop = (42 * storage->dim() + storage->dim() * storage->dim()) * storage->size() * storage->size();

	double needed_time = stopWatch.stop();
	double gflops = (flop / needed_time) / 1000000000;
	double gops = (gop / needed_time) / 1000000000;

	for(int i = 0; i < process_count ; i++)
	{
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if(i == process_index)
		{            
			//if (process_index == 0)
			std::cout << "[PROCESS :" << process_index << "] GFLOPS :" << gflops << " = (" << (flop/1000000000) << " / "<< needed_time << ")  GOPS :" << gops << std::endl;
		}
		std::cout.flush();
	}
#endif
}

}

}
