/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifdef USEMIC
#include "parallel/datadriven/basis/common/mic/MICKernel.hpp"
#include "parallel/datadriven/basis/common/mic/MICCPUHybridKernel.hpp"
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/MICLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/MICModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/MICModLinearMask.hpp"
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
#endif

#include "DMSystemMatrixMPITypeFactory.hpp"

#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinear.hpp"

#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMask.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAsyncMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityTrueAsyncMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityOnesidedMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAllreduce.hpp"

#include "parallel/datadriven/basis/common/CPUKernel.hpp"
#ifdef USEOCL
#include "parallel/datadriven/basis/common/ocl/OCLKernel.hpp"
#include "parallel/datadriven/basis/common/ocl/OCLKernelImpl.hpp"
#include "parallel/datadriven/basis/common/ocl/OCLCPUHybridKernel.hpp"

#include "parallel/datadriven/basis/linear/noboundary/operation/impl/OCLLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/OCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/OCLModLinearMask.hpp"
#endif

#include <cstring>
#include "base/exception/factory_exception.hpp"

#include "parallel/tools/MPI/SGppMPITools.hpp"

namespace sg {
  namespace parallel {
    template<typename KernelImplementation>
    datadriven::DMSystemMatrixBase* DMSystemMatrixMPITypeFactory::createDMSystemMatrixMPIType(base::Grid& grid, base::DataMatrix& trainDataset, double lambda, VectorizationType vecType, MPIType mpiType) {
      std::string parallelizationType;
      datadriven::DMSystemMatrixBase* result = 0;

      switch (mpiType) {
        case MPIAllreduce:
          parallelizationType = "Allreduce";
          result = new sg::parallel::DMSystemMatrixVectorizedIdentityAllreduce<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPIAlltoallv:
          parallelizationType = "Alltoallv";
          result = new sg::parallel::DMSystemMatrixVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
          break;

        case MPIAsync:
          parallelizationType = "Asynchronous Communication";
          result = new sg::parallel::DMSystemMatrixVectorizedIdentityAsyncMPI<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPITrueAsync:
          parallelizationType = "True Asynchronous Communication";
          result = new sg::parallel::DMSystemMatrixVectorizedIdentityTrueAsyncMPI<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPIOnesided:
          parallelizationType = "Onesided Communication";
          result = new sg::parallel::DMSystemMatrixVectorizedIdentityOnesidedMPI<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPINone:
          parallelizationType = "No MPI Implementation is used.";
          result = new sg::parallel::DMSystemMatrixVectorizedIdentity(grid, trainDataset, lambda, vecType);
          break;

        default:
          throw sg::base::factory_exception("this type of MPI communication is not yet implemented");
          break;
      }

      if (sg::parallel::myGlobalMPIComm->getMyRank() == 0) {
        std::cout << "Using MPI Parallelization: " << parallelizationType << " (" << sg::parallel::myGlobalMPIComm->getNumRanks() << " Processes)" << std::endl;
        size_t thread_count = 1;
#ifdef _OPENMP
        #pragma omp parallel
        {
          thread_count = omp_get_num_threads();
        }
#endif
        std::cout << "OpenMP: " << thread_count << " Threads active" << std::endl;
      }

      return result;
    }

    datadriven::DMSystemMatrixBase* DMSystemMatrixMPITypeFactory::getDMSystemMatrix(base::Grid& grid, base::DataMatrix& trainDataset, double lambda, VectorizationType vecType, MPIType mpiType) {
      const char* modlinear_mode = getenv("SGPP_MODLINEAR_EVAL");

      if (modlinear_mode == NULL) {
        modlinear_mode = "mask";
      }

      if (strcmp(grid.getType(), "linear") == 0 || strcmp(grid.getType(), "linearBoundary") == 0
          || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return createDMSystemMatrixMPIType<CPUKernel<X86SimdLinear> >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return createDMSystemMatrixMPIType<OCLKernel<OCLLinear<double> > >
                 (grid, trainDataset, lambda, vecType, mpiType);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return createDMSystemMatrixMPIType<OCLCPUHybridKernel<X86SimdLinear, OCLLinear<double > > >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return createDMSystemMatrixMPIType<MICKernel<MICLinear> >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return createDMSystemMatrixMPIType<MICCPUHybridKernel<X86SimdLinear, MICLinear> >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#endif
#endif
        else {
          std::cout << "WARNING: vectorization not implemented for this type of MPI-Communication. Using default (alltoallv). Please fix this in MPITypeFactory!!!!!! (" << __LINE__ << "in" << __FILE__ << ")" << std::endl;
          return new sg::parallel::DMSystemMatrixVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
        }
      } else if (strcmp(grid.getType(), "modlinear") == 0) {
        if (vecType == parallel::X86SIMD) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPIType<CPUKernel<X86SimdModLinear> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPIType<CPUKernel<X86SimdModLinearMask> >
                (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPIType<OCLKernel<OCLModLinear<double > > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPIType<OCLKernel<OCLModLinearMask<double> > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPIType<OCLCPUHybridKernel<X86SimdModLinear, OCLModLinear<double > > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPIType<OCLCPUHybridKernel<X86SimdModLinearMask, OCLModLinearMask<double > > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#endif

#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPIType<MICKernel<MICModLinear> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPIType<MICKernel<MICModLinearMask> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }
#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPIType<MICCPUHybridKernel<X86SimdModLinear, MICModLinear> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPIType<MICCPUHybridKernel<X86SimdModLinearMask, MICModLinearMask> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }
#endif
#endif
        else {
          std::cout << "WARNING: vectorization not implemented for this type of MPI-Communication. Using default (alltoallv). Please fix this in MPITypeFactory." << std::endl;
          return new sg::parallel::DMSystemMatrixVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
        }
      } else {
        throw base::factory_exception("MPITypeFactory: OperationMultipleEvalVectorized is not implemented for this grid type.");
      }
    }
  }
}
