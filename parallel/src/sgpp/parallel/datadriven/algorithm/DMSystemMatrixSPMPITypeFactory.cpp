// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifdef USEMIC
#include <sgpp/parallel/datadriven/basis/common/mic/SPMICKernel.hpp>
#include <sgpp/parallel/datadriven/basis/common/mic/SPMICCPUHybridKernel.hpp>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/impl/SPMICLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/SPMICModLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/SPMICModLinearMask.hpp>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
#endif

#include "DMSystemMatrixSPMPITypeFactory.hpp"

#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinear.hpp>

#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMask.hpp>

#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPI.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityAsyncMPI.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityTrueAsyncMPI.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityOnesidedMPI.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityAllreduce.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityBigdataAllreduce.hpp>

#include <sgpp/parallel/datadriven/basis/common/SPCPUKernel.hpp>
#ifdef USEOCL
#include <sgpp/parallel/datadriven/basis/common/ocl/SPOCLKernel.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/SPOCLKernelImpl.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/SPOCLCPUHybridKernel.hpp>

#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/impl/OCLLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/OCLModLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/OCLModLinearMask.hpp>
#endif

#include <cstring>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    template<typename KernelImplementation>
    datadriven::DMSystemMatrixBaseSP* DMSystemMatrixSPMPITypeFactory::createDMSystemMatrixMPITypeSP(base::Grid& grid, base::DataMatrixSP& trainDataset, float lambda, VectorizationType vecType, MPIType mpiType) {
      std::string parallelizationType;
      datadriven::DMSystemMatrixBaseSP* result = 0;

      switch (mpiType) {
        case MPIAllreduce:
          parallelizationType = "Allreduce";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityAllreduce<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPIAlltoallv:
          parallelizationType = "Alltoallv";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
          break;

        case MPIAsync:
          parallelizationType = "Asynchronous Communication";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityAsyncMPI<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPITrueAsync:
          parallelizationType = "True Asynchronous Communication";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityTrueAsyncMPI<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPIOnesided:
          parallelizationType = "Onesided Communication";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityOnesidedMPI<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPIBigdata:
          parallelizationType = "MPI Bigdata";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityBigdataAllreduce<KernelImplementation>(
            grid, trainDataset, lambda, vecType);
          break;

        case MPINone:
          parallelizationType = "No MPI Implementation is used.";
          result = new SGPP::parallel::DMSystemMatrixSPVectorizedIdentity(grid, trainDataset, lambda, vecType);
          break;

        default:
          throw SGPP::base::factory_exception("this type of MPI communication is not yet implemented");
          break;
      }

      if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
        std::cout << "Using MPI Parallelization: " << parallelizationType << " (" << SGPP::parallel::myGlobalMPIComm->getNumRanks() << " Processes)" << std::endl;
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

    datadriven::DMSystemMatrixBaseSP* DMSystemMatrixSPMPITypeFactory::getDMSystemMatrixSP(base::Grid& grid, base::DataMatrixSP& trainDataset, float lambda, VectorizationType vecType, MPIType mpiType) {
      const char* modlinear_mode = getenv("SGPP_MODLINEAR_EVAL");

      if (modlinear_mode == NULL) {
        modlinear_mode = "mask";
      }

      if (strcmp(grid.getType(), "linear") == 0 || strcmp(grid.getType(), "linearBoundary") == 0
          || strcmp(grid.getType(), "linearTruncatedBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return createDMSystemMatrixMPITypeSP<SPCPUKernel<SPX86SimdLinear> >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return createDMSystemMatrixMPITypeSP<SPOCLKernel<OCLLinear<float> > >
                 (grid, trainDataset, lambda, vecType, mpiType);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return createDMSystemMatrixMPITypeSP<SPOCLCPUHybridKernel<SPX86SimdLinear, OCLLinear<float > > >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return createDMSystemMatrixMPITypeSP<SPMICKernel<SPMICLinear> >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return createDMSystemMatrixMPITypeSP<SPMICCPUHybridKernel<SPX86SimdLinear, SPMICLinear> >
                 (grid, trainDataset, lambda, vecType, mpiType);
        }

#endif
#endif
        else {
          std::cout << "WARNING: vectorization not implemented for this type of MPI-Communication. Using default (alltoallv). Please fix this in MPITypeFactory!!!!!! (" << __LINE__ << "in" << __FILE__ << ")" << std::endl;
          return new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
        }
      } else if (strcmp(grid.getType(), "modlinear") == 0) {
        if (vecType == parallel::X86SIMD) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPITypeSP<SPCPUKernel<SPX86SimdModLinear> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPITypeSP<SPCPUKernel<SPX86SimdModLinearMask> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPITypeSP<SPOCLKernel<OCLModLinear<float > > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPITypeSP<SPOCLKernel<OCLModLinearMask<float> > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPITypeSP<SPOCLCPUHybridKernel<SPX86SimdModLinear, OCLModLinear<float > > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPITypeSP<SPOCLCPUHybridKernel<SPX86SimdModLinearMask, OCLModLinearMask<float > > >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#endif

#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPITypeSP<SPMICKernel<SPMICModLinear> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPITypeSP<SPMICKernel<SPMICModLinearMask> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return createDMSystemMatrixMPITypeSP<SPMICCPUHybridKernel<SPX86SimdModLinear, SPMICModLinear> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return createDMSystemMatrixMPITypeSP<SPMICCPUHybridKernel<SPX86SimdModLinearMask, SPMICModLinearMask> >
                   (grid, trainDataset, lambda, vecType, mpiType);
          } else {
            throw base::factory_exception("MPITypeFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#endif
#endif
        else {
          std::cout << "WARNING: vectorization not implemented for this type of MPI-Communication. Using default (alltoallv). Please fix this in MPITypeFactory." << std::endl;
          return new SGPP::parallel::DMSystemMatrixSPVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
        }
      } else {
        throw base::factory_exception("MPITypeFactory: OperationMultipleEvalVectorized is not implemented for this grid type.");
      }
    }
  }
}