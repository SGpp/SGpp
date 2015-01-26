#include "OCLPDEBound.hpp"

#ifdef USE_MPI
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#endif

namespace sg {
  namespace parallel {
    namespace oclpdekernels {

      cl_kernel ReduceBoundKernel[NUMDEVS];
      cl_mem d_ptrLevelBound[NUMDEVS];
      cl_mem d_ptrIndexBound[NUMDEVS];
      cl_mem d_ptrLevel_intBound[NUMDEVS];
      cl_mem d_ptrResultBound[NUMDEVS];
      cl_mem d_ptrResultPinnedBound[NUMDEVS];

      cl_mem d_ptrParResultBound[NUMDEVS];
      cl_mem d_ptrAlphaBound[NUMDEVS];
      cl_mem d_ptrLevelIndexLevelintconBound[NUMDEVS]; // constant memory buffer holding all three components
      cl_mem d_ptrLcl_qBound[NUMDEVS]; // Also holds q_inverse

      unsigned* offsetBound;
      REAL* ptrResultBound;

      REAL* ptrLevelTBound;
      REAL* ptrIndexTBound;
      REAL* ptrLevel_intTBound;
      REAL* ptrParResultBound;
      REAL* ptrAlphaEndBound;
      REAL* ptrLevelIndexLevelintBound;  // for the constant memory buffer holding all three components
      REAL* ptrLcl_qBound;             // Also holds q_inverse
      REAL* ptrResultTempBound;
      REAL* ptrResultZeroBound;
      REAL* ptrResultPinnedBound;


      size_t storageSizeBound;
      size_t storageSizePaddedBound;
      size_t num_groupsBound;
      size_t max_buffer_sizeBound;
      size_t par_result_sizeBound;
      size_t lcl_q_sizeBound;
      size_t alphaend_sizeBound;

      size_t storageInnerSizePaddedBound;
      size_t storageInnerSizeBound;
      size_t InnerSizeBound;
      size_t Inner_num_groupsBound;
      size_t Inner_par_result_sizeBound;
      size_t Inner_par_result_max_sizeBound;
      size_t Inner_result_sizeBound;
      size_t par_result_max_sizeBound;
      size_t pinnedResultBound_size;

      size_t constant_mem_size;
      size_t constant_buffer_size;
      size_t constant_buffer_iterations;

      size_t isFirstTimeLaplaceBound = 1;
      size_t isFirstTimeLTwoDotBound = 1;
      size_t isFirstTimeLTwoDotLaplaceBound = 1;
      // Used by MPI
#ifdef USE_MPI

      int* MPIOffsetListBound;
      int* MPISizeListBound;
#endif
      void boundarytransposer(REAL* sink, REAL* source, size_t dim1, size_t dim2, sg::base::GridStorage* storage) {

        unsigned k = 0;

        for (unsigned i = 0; i < dim2; i++) {
          sg::base::GridIndex* curPoint = (*storage)[i];

          if (curPoint->isInnerPoint()) {
            for (unsigned j = 0; j < dim1; j++) {
              sink[j * storageInnerSizePaddedBound + k] = source[i * dim1 + j];
            }

            k++;
          }
        }
      }


      std::string ReduceBoundKernelStr() {
        std::stringstream stream_program_src;
        stream_program_src << "__kernel void ReduceBoundKernel(				"  << std::endl;
        stream_program_src <<   "		   __global  REAL* ptrResultBound,	"  << std::endl;
        stream_program_src <<   "		   __global  REAL* ptrParResultBound,	"  << std::endl;
        stream_program_src <<   "		   ulong overallMultOffset,		"  << std::endl;
        stream_program_src <<   "		   ulong num_groups		"  << std::endl;
        stream_program_src <<   "		   )					"  << std::endl;
        stream_program_src <<   "{							"  << std::endl;
        stream_program_src <<   "unsigned j = get_global_id(0);				"  << std::endl;
        stream_program_src <<   "REAL res = 0.0;					"  << std::endl;
        stream_program_src <<   "for (unsigned k = 0; k < num_groups; k++) {		"  << std::endl;
        stream_program_src <<   "  res += ptrParResultBound[k*get_global_size(0) + j];	"  << std::endl;
        stream_program_src <<   "}							"  << std::endl;
        stream_program_src <<   "ptrResultBound[j] += res;		"  << std::endl;
        stream_program_src <<   "}							"  << std::endl;
        return stream_program_src.str();
      }

      void CompileReduceBound(int id, std::string kernel_src, cl_kernel* kernel) {

        cl_int err = CL_SUCCESS;
        std::string source2 = ReduceBoundKernelStr();

        std::stringstream stream_program_src;

        stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
        stream_program_src << "#define REAL double" << std::endl;
        stream_program_src << "#define LSIZE " << LSIZE << std::endl;
        stream_program_src << source2 << std::endl;
        std::string program_src = stream_program_src.str();
        const char* source3 = program_src.c_str();
        //  std::cout << "SOURCE " << source2 << std::endl;
        cl_program program = clCreateProgramWithSource(context, 1, (const char**)&source3, NULL, &err);
        oclCheckErr(err, "clCreateProgramWithSource");
        char buildOptions[256];
        int ierr = snprintf(buildOptions, sizeof(buildOptions), "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DINNERSTORAGEPAD=%lu -DDIMS=%lu -cl-finite-math-only -cl-fast-relaxed-math ", storageSizeBound, storageSizePaddedBound, num_groupsBound, storageInnerSizePaddedBound, dims);

        if (ierr < 0) {
          printf("Error in Build Options");
          exit(-1);
        }

        err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: compileReduceBound OpenCL Build Error. Error Code: " << err << std::endl;

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
      } // compileReduceBound


      std::string BoundLTwoDotFunction() {
        std::stringstream stream_program_src;
        stream_program_src << "REAL l2dot(	REAL lid, 	" << std::endl;
        stream_program_src << "   		REAL iid, 	" << std::endl;
        stream_program_src << "   		REAL in_lid, 	" << std::endl;
        stream_program_src << "   		REAL ljd, 	" << std::endl;
        stream_program_src << "   		REAL ijd,	" << std::endl;
        stream_program_src << "   		REAL in_ljd,	" << std::endl;
        stream_program_src << "   		REAL lcl_q)	" << std::endl;
        stream_program_src << "{				" << std::endl;
        stream_program_src << "  double res_one = select(0.0,(2.0/3.0) * in_lid, (ulong)((iid == ijd) && (ljd != 1)));	" << std::endl;
        stream_program_src << "  ulong selector = (lid > ljd);								" << std::endl;
        stream_program_src << "  double i1d = select(ijd, iid, selector);						" << std::endl;
        stream_program_src << "  double in_l1d = select(in_ljd, in_lid,selector);					" << std::endl;
        stream_program_src << "  double i2d = select(iid, ijd, selector);						" << std::endl;
        stream_program_src << "  double l2d = select(lid, ljd, selector);						" << std::endl;
        stream_program_src << "  double in_l2d = select(in_lid, in_ljd, selector);					" << std::endl;

#ifdef USEOCL_CPU
        stream_program_src << "  double q = (i1d * in_l1d -in_l1d); 							" << std::endl;
        stream_program_src << "  double p = (i1d * in_l1d +  in_l1d); 						     	" << std::endl;
        stream_program_src << "  ulong overlap = (max(q, (i2d * in_l2d - in_l2d)) < min(p, (i2d * in_l2d + in_l2d)));	" << std::endl;
        stream_program_src << "  double temp_res_inner = 2.0 - fabs((l2d * q - i2d)) - fabs((l2d * p - i2d));		" << std::endl;
#else
        stream_program_src << "  double q = fma(i1d, in_l1d, -in_l1d); 							" << std::endl;
        stream_program_src << "  double p = fma(i1d, in_l1d,  in_l1d); 						     	" << std::endl;
        stream_program_src << "  ulong overlap = (max(q, fma(i2d, in_l2d, -in_l2d)) < min(p, fma(i2d, in_l2d,in_l2d)));	" << std::endl;
        stream_program_src << "  double temp_res_inner = 2.0 - fabs(fma(l2d,q,-i2d)) - fabs(fma(l2d,p,-i2d));		" << std::endl;

#endif

        stream_program_src << "  double temp_res_rightbound = p + q;							" << std::endl;
        stream_program_src << "  double temp_res_leftbound = 2.0 - temp_res_rightbound;					" << std::endl;
        stream_program_src << "  double temp_res = select((temp_res_inner), 						" << std::endl;
        stream_program_src << "     select(temp_res_rightbound, temp_res_leftbound, (ulong)(i2d == 0))			" << std::endl;
        stream_program_src << "			   , (ulong) (l2d == 1));						" << std::endl;
        stream_program_src << "  temp_res *= (0.5 * in_l1d);								" << std::endl;
        stream_program_src << "  double res_two = select(0.0,temp_res,  overlap);					" << std::endl;
        stream_program_src << "  return (select(res_two, res_one, (ulong)(lid == ljd))) * lcl_q;			" << std::endl;
        stream_program_src << "}											" << std::endl;
        return stream_program_src.str();
      }

      void SetBuffersBound(REAL* ptrLevel,
                           REAL* ptrIndex,
                           REAL* ptrLevel_int,
                           size_t localStorageSize,
                           size_t localdim, sg::base::GridStorage* storage) {
        padding_size = num_devices * LSIZE;
        storageSizeBound = localStorageSize;
        size_t pad = padding_size - (storageSizeBound % padding_size);
        storageSizePaddedBound = storageSizeBound + pad;
        dims = localdim;
        cl_ulong sizec;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &sizec, 0);
        constant_mem_size = sizec;

        num_groupsBound = (storageSizePaddedBound) / LSIZE / num_devices;
        size_t level_size = storageSizePaddedBound * dims;
        size_t index_size = storageSizePaddedBound * dims;
        size_t level_int_size = storageSizePaddedBound * dims;
        size_t alpha_size = storageSizePaddedBound;
        size_t lambda_size = dims;
        par_result_sizeBound = storageSizePaddedBound * num_groupsBound * num_devices;
        lcl_q_sizeBound      = dims + dims;
        constant_mem_size = constant_mem_size - lcl_q_sizeBound * sizeof(REAL)
                            - lambda_size * sizeof(REAL);
        alphaend_sizeBound = pad;


#if PRINTBUFFERSIZES
        std::cout << "storageSizeBound " << storageSizeBound << std::endl;
        std::cout << "storageSizePaddedBound " << storageSizePaddedBound << std::endl;
        std::cout << "dims " << dims << std::endl;
        std::cout << "num_groupsBound " << num_groupsBound << std::endl;
        std::cout << "index_size " << index_size << std::endl;
        std::cout << "level_int_size " << level_int_size << std::endl;
        std::cout << "alpha_size " << alpha_size << std::endl;
        std::cout << "par_result_sizeBound " << par_result_sizeBound << std::endl;
        std::cout << "lcl_q_sizeBound " << lcl_q_sizeBound << std::endl;
#endif
        //       std::cout << level_size << " " << index_size << std::endl;
        cl_int ciErrNum = CL_SUCCESS;

        ptrLcl_qBound = (REAL*)calloc(lcl_q_sizeBound, sizeof(REAL));
        ptrAlphaEndBound = (REAL*)calloc(alphaend_sizeBound, sizeof(REAL));

        size_t innerpoints = storage->getNumInnerPoints();
        size_t pad2 = padding_size - (innerpoints % padding_size);
        storageInnerSizeBound = innerpoints;
        storageInnerSizePaddedBound = innerpoints + pad2;

        InnerSizeBound = storageInnerSizePaddedBound * dims;
        Inner_num_groupsBound = (storageSizePaddedBound) / LSIZE / num_devices;
        Inner_par_result_sizeBound = storageInnerSizePaddedBound * Inner_num_groupsBound * num_devices;
        Inner_result_sizeBound = storageInnerSizePaddedBound;
        pinnedResultBound_size = storageInnerSizePaddedBound * num_devices;


        ptrLevelTBound = (REAL*)calloc(InnerSizeBound, sizeof(REAL));
        ptrIndexTBound = (REAL*)calloc(InnerSizeBound, sizeof(REAL));
        ptrLevel_intTBound = (REAL*)calloc(InnerSizeBound, sizeof(REAL));
        offsetBound = (unsigned*)calloc(storageInnerSizePaddedBound, sizeof(unsigned));
        ptrResultBound = (REAL*)calloc(Inner_result_sizeBound, sizeof(REAL));
        ptrResultTempBound = (REAL*)calloc(pinnedResultBound_size, sizeof(REAL));
        ptrResultZeroBound = (REAL*)calloc(Inner_result_sizeBound, sizeof(REAL));

        unsigned num = 0;

        for (unsigned i = 0; i < storage->size(); i++) {
          sg::base::GridIndex* curPoint = (*storage)[i];

          if (curPoint->isInnerPoint()) {
            offsetBound[num] = i;
            num++;
          }
        }


        boundarytransposer(ptrLevelTBound, ptrLevel, dims, storageSizeBound, storage);
        boundarytransposer(ptrIndexTBound, ptrIndex, dims, storageSizeBound, storage);
        boundarytransposer(ptrLevel_intTBound, ptrLevel_int, dims, storageSizeBound, storage);


        cl_ulong size3;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size3, 0);

        cl_ulong size4;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_MAX_MEM_ALLOC_SIZE , sizeof(cl_ulong), &size4, 0);

        size_t sizedoubles = size3 / sizeof(REAL);
        size_t sizedoubles64 = sizedoubles;
        size_t gpu_max_buffer_size = size4 / sizeof(REAL);
        gpu_max_buffer_size = gpu_max_buffer_size - (gpu_max_buffer_size % (storageSizePaddedBound));

        size_t memoryOther = Inner_result_sizeBound + 3 * InnerSizeBound + lcl_q_sizeBound;
        size_t memoryNonParResult = level_size + index_size + level_int_size + alpha_size + lambda_size;
        size_t memoryLeftover = sizedoubles64 - memoryOther - memoryNonParResult;
        Inner_par_result_max_sizeBound = memoryLeftover - (memoryLeftover % storageSizePaddedBound);

        max_buffer_sizeBound = MAXBYTES / sizeof(REAL) - memoryOther - memoryNonParResult;
        max_buffer_sizeBound = max_buffer_sizeBound - ((max_buffer_sizeBound) % (storageSizePaddedBound));

        //  std::cout << " max_buffer_sizeBound " << max_buffer_sizeBound <<std::endl;
        Inner_par_result_max_sizeBound = std::min(Inner_par_result_max_sizeBound, max_buffer_sizeBound);

        Inner_par_result_max_sizeBound = std::min(Inner_par_result_max_sizeBound, Inner_par_result_sizeBound);
        Inner_par_result_max_sizeBound = std::min(Inner_par_result_max_sizeBound, gpu_max_buffer_size) / num_devices;
        //  std::cout << " gpu_max_buffer_size " << gpu_max_buffer_size<<std::endl;

        ptrParResultBound = (REAL*)calloc(Inner_par_result_max_sizeBound, sizeof(REAL));

        ptrLevelIndexLevelintBound = (REAL*)calloc(level_size + index_size + level_int_size, sizeof(REAL));
        size_t three_d = 0;
        size_t threedims = 3 * dims;

        for (size_t i = 0; i < storageSizeBound; i++) {
          for (size_t d = 0; d < dims; d++) {
            ptrLevelIndexLevelintBound[i * threedims + three_d] = ptrLevel[i * dims + d];
            three_d += 1;
            ptrLevelIndexLevelintBound[i * threedims + three_d] = ptrIndex[i * dims + d];
            three_d += 1;
            ptrLevelIndexLevelintBound[i * threedims + three_d] = ptrLevel_int[i * dims + d];
            three_d += 1;
          }

          three_d = 0;
        }

        constant_buffer_size = (constant_mem_size / sizeof(REAL));
        size_t mod_res = constant_buffer_size % (3 * dims * padding_size);
        constant_buffer_size = constant_buffer_size - mod_res;
        constant_buffer_iterations = constant_buffer_size / (3 * dims);

        for (size_t i = 0; i < num_devices; ++i) {
          d_ptrLevelBound[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              InnerSizeBound * sizeof(REAL), ptrLevelTBound, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrLevel");
          d_ptrLevelIndexLevelintconBound[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                               constant_buffer_size * sizeof(REAL), NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrLevelIndexLevelint");

          d_ptrIndexBound[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              InnerSizeBound * sizeof(REAL), ptrIndexTBound, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrIndex");

          d_ptrLevel_intBound[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                                  InnerSizeBound * sizeof(REAL), ptrLevel_intTBound, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer ptrLevel_int");

          d_ptrAlphaBound[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                              alpha_size * sizeof(REAL), NULL, &ciErrNum);


          oclCheckErr(ciErrNum, "clCreateBuffer ptrAlpha");




          //RED
          d_ptrResultPinnedBound[i] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR ,
                                      pinnedResultBound_size * sizeof(REAL), NULL, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer d_ptrResultPinned");

          d_ptrResultBound[i] = clCreateBuffer(context, CL_MEM_READ_WRITE ,
                                               Inner_result_sizeBound * sizeof(REAL), NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clCreateBuffer d_ptrResultBound");

          ptrResultPinnedBound = (REAL*) clEnqueueMapBuffer(command_queue[0],
                                 d_ptrResultPinnedBound[0],
                                 CL_TRUE,
                                 CL_MAP_READ,
                                 0,
                                 pinnedResultBound_size * sizeof(REAL),
                                 0, NULL, NULL, &ciErrNum);
          oclCheckErr(ciErrNum, "clEnqueueMapBuffer ptrResultPinnedBound");



          d_ptrParResultBound[i] = clCreateBuffer(context, CL_MEM_READ_WRITE ,
                                                  Inner_par_result_max_sizeBound * sizeof(REAL), NULL, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrParResult");

          d_ptrLcl_qBound[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                              lcl_q_sizeBound * sizeof(REAL), NULL, &ciErrNum);

          oclCheckErr(ciErrNum, "clCreateBuffer ptrLcl_q");

        }

      } // SetBuffersBound

      void CleanUpBound() {
        if (!isFirstTimeLaplaceBound || !isFirstTimeLTwoDotBound
            || !isFirstTimeLTwoDotLaplaceBound) {
          cl_int ciErrNum = CL_SUCCESS;

          for (unsigned int i = 0; i < num_devices; i++) {
            ciErrNum |= clReleaseMemObject(d_ptrLevelBound[i]);
            ciErrNum |= clReleaseMemObject(d_ptrIndexBound[i]);
            ciErrNum |= clReleaseMemObject(d_ptrLevel_intBound[i]);
            ciErrNum |= clReleaseMemObject(d_ptrResultBound[i]);
            ciErrNum |= clReleaseMemObject(d_ptrParResultBound[i]);
            ciErrNum |= clReleaseMemObject(d_ptrLevelIndexLevelintconBound[i]);

            ciErrNum |= clReleaseMemObject(d_ptrAlphaBound[i]);
            ciErrNum |= clReleaseMemObject(d_ptrLcl_qBound[i]);
            oclCheckErr(ciErrNum, "clReleaseMemObject");
            ciErrNum |= clReleaseKernel( ReduceBoundKernel[i] );
            oclCheckErr(ciErrNum, "clReleaseKernel");
          }

          free(ptrLevelTBound);
          free(ptrIndexTBound);
          free(ptrLevel_intTBound);
          free(ptrParResultBound);
          free(ptrLcl_qBound);
          free(ptrResultBound);
          free(offsetBound);


          CleanUpLaplaceBound();
          CleanUpLTwoDotBound();
          CleanUpLTwoDotLaplaceBound();
        }
      } // CleanUpBound

#ifdef USE_MPI
#ifdef USEOCL
      void SetUpMPIBound() {
        int myrank, nproz2;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproz2);

        size_t nproz = static_cast<size_t>(nproz2);
        MPIOffsetListBound = new int[nproz];
        MPISizeListBound = new int[nproz];

        size_t minimum_size  = num_devices * LSIZE;
        size_t data_size = storageSizePaddedBound;

        sg::parallel::PartitioningTool::calcDistribution(data_size, nproz, MPISizeListBound, MPIOffsetListBound, minimum_size);

        if (myrank == 0) {
          std::cout << "nproz " << nproz << std::endl;

          for (size_t i = 0; i < nproz; i++) {
            std::cout << "MPIOffsetListBound[" << i << "] = " << MPIOffsetListBound[i] << std::endl;
            std::cout << "MPISizeListBound[" << i << "] = " << MPISizeListBound[i] << std::endl;
          }

          std::cout << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);


      }

      void MPI_ShareResultAllReduceBound(sg::base::DataVector& result) {
        REAL* ptrResult = result.getPointer();
        MPI_Allreduce(MPI_IN_PLACE , ptrResultBound,
                      (int)(storageInnerSizeBound),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (size_t i = 0; i < storageInnerSizeBound ; i++) {
          ptrResult[offsetBound[i]] = ptrResultBound[i];
        }

      }
#endif
#endif

    }
  }
}
