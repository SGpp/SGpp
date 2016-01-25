// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OCLPDEInner.hpp"

#ifdef USE_MPI
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/tools/MPI/MPICommunicator.hpp>
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {
      cl_kernel ReduceInnerKernel[NUMDEVS];
      cl_mem d_ptrLevelInner[NUMDEVS];
      cl_mem d_ptrIndexInner[NUMDEVS];
      cl_mem d_ptrLevel_intInner[NUMDEVS];
      cl_mem d_ptrResultInner[NUMDEVS];
      cl_mem d_ptrResultPinnedInner[NUMDEVS];


      cl_mem d_ptrParResultInner[NUMDEVS];
      cl_mem d_ptrAlphaInner[NUMDEVS];
      cl_mem d_ptrAlphaPinnedInner[NUMDEVS];
      cl_mem d_ptrLevelIndexLevelintconInner[NUMDEVS]; // constant memory buffer holding all three components
      cl_mem d_ptrLcl_qInner[NUMDEVS]; // Also holds q_inverse

      REAL* ptrLevelTInner;
      REAL* ptrIndexTInner;
      REAL* ptrLevel_intTInner;
      REAL* ptrParResultInner;
      REAL* ptrAlphaEndInner;
      REAL* ptrLevelIndexLevelintInner;  // for the constant memory buffer holding all three components
      REAL* ptrLcl_qInner;             // Also holds q_inverse
      REAL* ptrResultTemp;
      REAL* ptrResultZero;
      REAL* ptrResultPinnedInner;
      REAL* ptrAlphaPinnedInner;

      size_t storageSize;
      size_t storageSizePadded;
      size_t num_groups;
      size_t par_result_size;
      size_t lcl_q_size;
      size_t alphaend_size;
      size_t max_buffer_size;
      size_t par_result_max_size;

      size_t constant_mem_size_noboundary;
      size_t constant_buffer_size_noboundary;
      size_t constant_buffer_iterations_noboundary;

      size_t isFirstTimeLaplaceInner = 1;
      size_t isFirstTimeLTwoDotInner = 1;
      size_t isFirstTimeLTwoDotLaplaceInner = 1;



      // Used by MPI
#ifdef USE_MPI

      int* MPIOffsetListInner;
      int* MPISizeListInner;
      MPICommunicator* MPICommunicator;
#endif



      void transposer(REAL* sink, REAL* source, size_t dim1, size_t dim2pad, size_t dim2) {
        for (unsigned i = 0; i < dim2; i++) {
          for (unsigned j = 0; j < dim1; j++) {
            sink[j * dim2pad + i] = source[i * dim1 + j];
          }
        }
      }

      std::string ReduceInnerKernelStr() {
        std::stringstream stream_program_src;
        stream_program_src <<   "__kernel void ReduceInnerKernel(__global  REAL* ptrResult,    "  << std::endl;
        stream_program_src << "                                __global  REAL* ptrParResult, "  << std::endl;
        stream_program_src << "                		 ulong overallParOffset,      "  << std::endl;
        stream_program_src << "                		 ulong num_groups)      "  << std::endl;
        stream_program_src << "{							       "  << std::endl;
        stream_program_src << "unsigned j = get_global_id(0);				       "  << std::endl;
        stream_program_src << "REAL res = 0.0;					       "  << std::endl;
        stream_program_src << "for (unsigned k = 0; k < num_groups; k++) {		       "  << std::endl;
        stream_program_src << "	res += ptrParResult[k*get_global_size(0) + j];         "  << std::endl;
        stream_program_src << "}							       "  << std::endl;
        stream_program_src << "ptrResult[j] += res;		       "  << std::endl;
        stream_program_src << "}							       "  << std::endl;
        return stream_program_src.str();
      }


      void CompileReduceInner(int id, std::string kernel_src, cl_kernel* kernel) {

        cl_int err = CL_SUCCESS;
        std::string source2 = ReduceInnerKernelStr();

        std::stringstream stream_program_src;

#ifdef USEOCL_NVIDIA
        stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
#endif
        stream_program_src << "#define REAL double" << std::endl;
        stream_program_src << "#define LSIZE " << LSIZE << std::endl;
        stream_program_src << source2 << std::endl;
        std::string program_src = stream_program_src.str();
        const char* source3 = program_src.c_str();
        //  std::cout << "SOURCE " << source2 << std::endl;
        cl_program program = clCreateProgramWithSource(context, 1, (const char**)&source3, NULL, &err);
        oclCheckErr(err, "clCreateProgramWithSource");
        char buildOptions[256];
        int ierr = snprintf(buildOptions, sizeof(buildOptions), "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DINNERSTORAGEPAD=%lu -DDIMS=%lu -cl-finite-math-only -cl-fast-relaxed-math ", storageSize, storageSizePadded, num_groups, 0L, dims);

        if (ierr < 0) {
          printf("Error in Build Options");
          exit(-1);
        }

        err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: compileReduceInner OpenCL Build Error. Error Code: " << err << std::endl;

          size_t len;
          char buffer[10000];

          // get the build log
          clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);

          std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
        }

        oclCheckErr(err, "clBuildProgram");


        kernel[id] = clCreateKernel(program, kernel_src.c_str(), &err);

        oclCheckErr(err, "clCreateKernel");

        err |= clReleaseProgram(program);
        oclCheckErr(err, "clReleaseProgram");
      }


      std::string InnerLTwoDotFunction() {
        std::stringstream stream_program_src;
        stream_program_src << "REAL l2dot(REAL lid, 	   	" << std::endl;
        stream_program_src << "	        REAL iid,		" << std::endl;
        stream_program_src << "		REAL in_lid,		" << std::endl;
        stream_program_src << "		REAL ljd,		" << std::endl;
        stream_program_src << "		REAL ijd,		" << std::endl;
        stream_program_src << "		REAL in_ljd,		" << std::endl;
        stream_program_src << "		REAL lcl_q)		" << std::endl;
        stream_program_src << "{					" << std::endl;
        stream_program_src << "double res_one = select(0.0,(2.0/3.0) * in_lid, (ulong)(iid == ijd));	" << std::endl;
        stream_program_src << "ulong selector = (lid > ljd);						" << std::endl;
        stream_program_src << "double i1d = select(ijd, iid, selector);				" << std::endl;
        stream_program_src << "double in_l1d = select(in_ljd, in_lid,selector);			" << std::endl;
        stream_program_src << "double i2d = select(iid, ijd, selector);				" << std::endl;
        stream_program_src << "double l2d = select(lid, ljd, selector);				" << std::endl;
        stream_program_src << "double in_l2d = select(in_lid, in_ljd, selector);			" << std::endl;
#ifdef USEOCL_CPU
        stream_program_src << "double q = (i1d-1)*in_l1d;		" << std::endl;
        stream_program_src << "double p = (i1d+1)*in_l1d;		" << std::endl;
        stream_program_src << "ulong overlap = (max(q, (i2d * in_l2d -in_l2d)) < min(p, (i2d * in_l2d + in_l2d)));			" << std::endl;
        stream_program_src << "double temp_res = ((0.5*in_l1d) * (- fabs((l2d *q -i2d)) - fabs((l2d * p - i2d))) + in_l1d);	" << std::endl;
#else
        stream_program_src << "double q = fma(i1d, in_l1d, -in_l1d); //(i1d-1)*in_l1d;		" << std::endl;
        stream_program_src << "double p = fma(i1d, in_l1d,  in_l1d); //(i1d+1)*in_l1d;		" << std::endl;
        stream_program_src << "ulong overlap = (max(q, fma(i2d, in_l2d, -in_l2d)) < min(p, fma(i2d, in_l2d,in_l2d)));			" << std::endl;
        stream_program_src << "double temp_res = fma((0.5*in_l1d), (- fabs(fma(l2d,q,-i2d)) - fabs(fma(l2d,p,-i2d))), in_l1d);	" << std::endl;
#endif
        stream_program_src << "double res_two = select(0.0,temp_res,overlap); // Now mask result					" << std::endl;
        stream_program_src << "return (select(res_two, res_one, (ulong)(lid == ljd)))*lcl_q;						" << std::endl;
        stream_program_src << "}"   << std::endl;

        return stream_program_src.str();
      }


      void SetBuffersInner(REAL* ptrLevel,
                           REAL* ptrIndex,
                           REAL* ptrLevel_int,
                           size_t localStorageSize,
                           size_t localdim, SGPP::base::GridStorage* storage) {
        padding_size = num_devices * LSIZE;

        storageSize = localStorageSize;
        size_t pad = padding_size - (storageSize % padding_size);
        storageSizePadded = storageSize + pad;
        dims = localdim;
        cl_ulong sizec;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &sizec, 0);
        constant_mem_size_noboundary = sizec;


        num_groups = (storageSizePadded) / LSIZE / num_devices;
        size_t level_size = storageSizePadded * dims;
        size_t index_size = storageSizePadded * dims;
        size_t level_int_size = storageSizePadded * dims;
        size_t alpha_size = storageSizePadded;
        size_t lambda_size = dims;
        size_t result_size = storageSizePadded;
        size_t pinnedresult_size = storageSizePadded * num_devices;
        par_result_size = storageSizePadded * num_groups * num_devices;
        lcl_q_size      = dims + dims;
        constant_mem_size_noboundary = constant_mem_size_noboundary - lcl_q_size * sizeof(REAL)
                                       - lambda_size * sizeof(REAL);
        alphaend_size = pad;
#ifdef USE_MPI
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

#ifdef USE_MPI

        if (myrank == 0) {
#endif
#if PRINTBUFFERSIZES
          std::cout << "storageSize " << storageSize << std::endl;
          std::cout << "storageSizePadded " << storageSizePadded << std::endl;
          std::cout << "dims " << dims << std::endl;
          std::cout << "num_groups " << num_groups << std::endl;
          std::cout << "level_size " << level_size << std::endl;
          std::cout << "index_size " << index_size << std::endl;
          std::cout << "level_int_size " << level_int_size << std::endl;
          std::cout << "alpha_size " << alpha_size << std::endl;
          std::cout << "result_size " << result_size << std::endl;
          std::cout << "par_result_size " << par_result_size << std::endl;
          std::cout << "lcl_q_size " << lcl_q_size << std::endl;
#endif
#ifdef USE_MPI
        }

#endif

        //       std::cout << level_size << " " << index_size << std::endl;
        cl_int ciErrNum = CL_SUCCESS;
        ptrLevelTInner = (REAL*)calloc(level_size, sizeof(REAL));
        ptrIndexTInner = (REAL*)calloc(index_size, sizeof(REAL));
        ptrLevel_intTInner = (REAL*)calloc(level_int_size, sizeof(REAL));
        transposer(ptrLevelTInner, ptrLevel, dims, storageSizePadded, storageSize);
        transposer(ptrIndexTInner, ptrIndex, dims, storageSizePadded, storageSize);
        transposer(ptrLevel_intTInner, ptrLevel_int, dims, storageSizePadded, storageSize);


        ptrResultTemp = (REAL*)calloc(pinnedresult_size, sizeof(REAL));
        ptrResultZero = (REAL*)calloc(result_size, sizeof(REAL));
        ptrLcl_qInner = (REAL*)calloc(lcl_q_size, sizeof(REAL));
        ptrAlphaEndInner = (REAL*)calloc(alphaend_size, sizeof(REAL));

        ptrLevelIndexLevelintInner = (REAL*)calloc(level_size + index_size + level_int_size, sizeof(REAL));

        size_t three_d = 0;
        size_t threedims = 3 * dims;

        for (size_t i = 0; i < storageSize; i++) {
          for (size_t d = 0; d < dims; d++) {
            ptrLevelIndexLevelintInner[i * threedims + three_d] = ptrLevel[i * dims + d];
            three_d += 1;
            ptrLevelIndexLevelintInner[i * threedims + three_d] = ptrIndex[i * dims + d];
            three_d += 1;
            ptrLevelIndexLevelintInner[i * threedims + three_d] = ptrLevel_int[i * dims + d];
            three_d += 1;
          }

          three_d = 0;
        }

        cl_ulong size3;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size3, 0);

        cl_ulong size4;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_MAX_MEM_ALLOC_SIZE , sizeof(cl_ulong), &size4, 0);


        size_t sizedoubles64 = size3 / sizeof(REAL);
        size_t gpu_max_buffer_size = size4 / sizeof(REAL);
        gpu_max_buffer_size = gpu_max_buffer_size - (gpu_max_buffer_size % storageSizePadded);
        size_t memoryNonParResult = level_size + index_size + level_int_size + alpha_size + result_size;
        size_t memoryLeftover = sizedoubles64 - memoryNonParResult;
        par_result_max_size = memoryLeftover - (memoryLeftover % storageSizePadded);
        par_result_max_size = std::min(par_result_max_size, par_result_size);
        par_result_max_size = std::min(par_result_max_size, gpu_max_buffer_size) / num_devices;



        // #ifdef USE_MPI
        //  if(myrank == 0) {
        //    std::cout << " par_result_max_size " << par_result_max_size << std::endl;
        //  }
        // #endif
        //      std::cout << " par_result_max_size " << par_result_max_size << std::endl;

        ptrParResultInner = (REAL*)calloc(par_result_max_size, sizeof(REAL));


        constant_buffer_size_noboundary = (constant_mem_size_noboundary / sizeof(REAL));

        size_t mod_res = constant_buffer_size_noboundary % (3 * dims * padding_size);
        constant_buffer_size_noboundary = constant_buffer_size_noboundary - mod_res;
        constant_buffer_iterations_noboundary = constant_buffer_size_noboundary / (3 * dims);


        for (size_t i = 0; i < num_devices; ++i) {

          d_ptrLevelInner[i] = clCreateBuffer(context, CL_MEM_READ_ONLY
                                              | CL_MEM_COPY_HOST_PTR
                                              ,
                                              level_size * sizeof(REAL), ptrLevelTInner, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrLevelInnerL290");


          d_ptrLevelIndexLevelintconInner[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                               constant_buffer_size_noboundary * sizeof(REAL), NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrLevelIndexLevelint");

          d_ptrIndexInner[i] = clCreateBuffer(context,
                                              CL_MEM_READ_ONLY
                                              | CL_MEM_COPY_HOST_PTR
                                              ,
                                              index_size * sizeof(REAL),
                                              ptrIndexTInner, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrIndex");


          d_ptrLevel_intInner[i] = clCreateBuffer(context, CL_MEM_READ_ONLY
                                                  | CL_MEM_COPY_HOST_PTR
                                                  ,
                                                  level_int_size * sizeof(REAL), ptrLevel_intTInner , &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrLevel_int");

          d_ptrAlphaInner[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                              alpha_size * sizeof(REAL), NULL, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrAlpha");

          d_ptrAlphaPinnedInner[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR ,
                                     alpha_size * sizeof(REAL), NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrAlphaPinnedInner");

          ptrAlphaPinnedInner = (REAL*) clEnqueueMapBuffer(command_queue[0],
                                d_ptrAlphaPinnedInner[0],
                                CL_TRUE,
                                CL_MAP_WRITE,
                                0,
                                alpha_size * sizeof(REAL),
                                0, NULL, NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clEnqueueMapBuffer ptrAlphaPinned");

          for (size_t jj = 0; jj < storageSizePadded; jj++) {
            ptrAlphaPinnedInner[jj] = 0.0;
          }

          //RED
          d_ptrResultPinnedInner[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR ,
                                      pinnedresult_size * sizeof(REAL), NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrResultPinnedInner");
          d_ptrResultInner[i] = clCreateBuffer(context, CL_MEM_READ_WRITE ,
                                               result_size * sizeof(REAL), NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrResult");

          ptrResultPinnedInner = (REAL*) clEnqueueMapBuffer(command_queue[0],
                                 d_ptrResultPinnedInner[0],
                                 CL_TRUE,
                                 CL_MAP_READ,
                                 0,
                                 pinnedresult_size * sizeof(REAL),
                                 0, NULL, NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clEnqueueMapBuffer ptrResult");



          d_ptrParResultInner[i] = clCreateBuffer(context, CL_MEM_READ_WRITE ,
                                                  par_result_max_size * sizeof(REAL), NULL, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrParResult");

          d_ptrLcl_qInner[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                              lcl_q_size * sizeof(REAL), NULL, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrLcl_q");

        }
      } // SetBuffersInner

      void CleanUpInner() {
        if (!isFirstTimeLaplaceInner || !isFirstTimeLTwoDotInner
            || !isFirstTimeLTwoDotLaplaceInner) {
          cl_int ciErrNum = CL_SUCCESS;

          for (unsigned int i = 0; i < num_devices; i++) {
            ciErrNum |= clReleaseMemObject(d_ptrLevelInner[i]);
            ciErrNum |= clReleaseMemObject(d_ptrIndexInner[i]);
            ciErrNum |= clReleaseMemObject(d_ptrLevel_intInner[i]);
            ciErrNum |= clReleaseMemObject(d_ptrResultInner[i]);
            ciErrNum |= clReleaseMemObject(d_ptrParResultInner[i]);
            ciErrNum |= clReleaseMemObject(d_ptrLevelIndexLevelintconInner[i]);

            ciErrNum |= clReleaseMemObject(d_ptrAlphaInner[i]);
            ciErrNum |= clReleaseMemObject(d_ptrLcl_qInner[i]);
            ciErrNum |= clReleaseKernel( ReduceInnerKernel[i] );
            oclCheckErr(ciErrNum, "clReleaseMemObject & clReleaseKernel");
          }

          free(ptrLevelTInner);
          free(ptrIndexTInner);
          free(ptrLevel_intTInner);
          free(ptrParResultInner);
          free(ptrLcl_qInner);

          CleanUpLaplaceInner();
          CleanUpLTwoDotInner();
          CleanUpLTwoDotLaplaceInner();
        }
      }

#ifdef USE_MPI
#ifdef USEOCL
      void SetUpMPIInner() {
        int myrank, nproz2;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproz2);

        size_t nproz = static_cast<size_t>(nproz2);

        MPIOffsetListInner = new int[nproz];
        MPISizeListInner = new int[nproz];

        size_t minimum_size  = num_devices * LSIZE;
        size_t data_size = storageSizePadded;

        SGPP::parallel::PartitioningTool::calcDistribution(data_size, nproz, MPISizeListInner, MPIOffsetListInner, minimum_size);

        if (myrank == 0) {
          std::cout << "nproz " << nproz << std::endl;

          for (size_t i = 0; i < nproz; i++) {
            std::cout << "MPIOffsetListInner[" << i << "] = " << MPIOffsetListInner[i] << std::endl;
            std::cout << "MPISizeListInner[" << i << "] = " << MPISizeListInner[i] << std::endl;
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPICommunicator = new SGPP::parallel::MPICommunicator(myrank, nproz2);

      }

      void MPI_CombineResultInner(SGPP::base::DataVector& result) {
        MPICommunicator->reduceGridCoefficients(result);
        //  double * ptrResult = result.getPointer();
        //  MPI_Allreduce(MPI_IN_PLACE, ptrResult,
        //          storageSize, MPI_DOUBLE,
        //          MPI_SUM, MPI_COMM_WORLD);
      }
#endif
#endif
    }
  }
}
