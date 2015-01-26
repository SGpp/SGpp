
#ifdef USE_MPI
#ifdef USEOCL
#include <mpi.h>
#endif
#endif
#include "OCLPDEKernels.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      extern cl_kernel ReduceInnerKernel[NUMDEVS];
      extern cl_mem d_ptrLevelInner[NUMDEVS];
      extern cl_mem d_ptrIndexInner[NUMDEVS];
      extern cl_mem d_ptrLevel_intInner[NUMDEVS];
      extern cl_mem d_ptrResultInner[NUMDEVS];
      extern cl_mem d_ptrResultPinnedInner[NUMDEVS];

      extern cl_mem d_ptrParResultInner[NUMDEVS];
      extern cl_mem d_ptrAlphaInner[NUMDEVS];
      extern cl_mem d_ptrAlphaPinnedInner[NUMDEVS];
      extern cl_mem d_ptrLevelIndexLevelintconInner[NUMDEVS]; // constant memory buffer holding all three components
      extern cl_mem d_ptrLcl_qInner[NUMDEVS]; // Also holds q_inverse

      extern REAL* ptrLevelTInner;
      extern REAL* ptrIndexTInner;
      extern REAL* ptrLevel_intTInner;
      extern REAL* ptrParResultInner;
      extern REAL* ptrAlphaEndInner;
      extern REAL* ptrLevelIndexLevelintInner;  // for the constant memory buffer holding all three components
      extern REAL* ptrLcl_qInner;             // Also holds q_inverse
      extern REAL* ptrResultTemp;
      extern REAL* ptrResultZero;
      extern REAL* ptrResultPinnedInner;
      extern REAL* ptrAlphaPinnedInner;


      extern size_t storageSize;
      extern size_t storageSizePadded;
      extern size_t num_groups;
      extern size_t par_result_size;
      extern size_t lcl_q_size;
      extern size_t alphaend_size;
      extern size_t max_buffer_size;
      extern size_t par_result_max_size;

      extern size_t constant_mem_size_noboundary;
      extern size_t constant_buffer_size_noboundary;
      extern size_t constant_buffer_iterations_noboundary;

      extern size_t isFirstTimeLaplaceInner;
      extern size_t isFirstTimeLTwoDotInner;
      extern size_t isFirstTimeLTwoDotLaplaceInner;

      /// Returns the string with the OpenCL code for the Reduction kernel for the operators on the inner grid.
      std::string ReduceInnerKernelStr();
      /// Compiles the OpenCL code for the Reduction kernel for the operators on the inner grid and saves it in kernel[id]. kernel_src must match the name of the OpenCL function.
      void CompileReduceInner(int id, std::string kernel_src, cl_kernel* kernel);

      /// Returns the string with the OpenCL code for the LTwoDot function for the operators on the inner grid.
      std::string InnerLTwoDotFunction();

      /// Allocates and initializes the main part of the buffers needed by the OpenCL code for the operators on the inner grid.
      void SetBuffersInner(REAL* ptrLevel,
                           REAL* ptrIndex,
                           REAL* ptrLevel_int,
                           size_t localStorageSize,
                           size_t localdim, SGPP::base::GridStorage* storage);
#ifdef USE_MPI
      extern int* MPIOffsetListInner;
      extern int* MPISizeListInner;
      /// Creates partitions of the work load along the j-direction for the operators on the inner grid.
      void SetUpMPIInner();
      /// Computes the sum of all subresults result on all processes.
      void MPI_CombineResultInner(SGPP::base::DataVector& result);
#endif
      /// Deallocates all data pertaining to the Laplace operator on the inner grid
      void CleanUpLaplaceInner();
      /// Deallocates all data pertaining to the LTwoDot operator on the inner grid
      void CleanUpLTwoDotInner();
      /// Deallocates all data pertaining to the combined LTwoDot+Laplace Operator working on the inner grid
      void CleanUpLTwoDotLaplaceInner();
    }
  }
}
