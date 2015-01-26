// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYMPIBASE_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYMPIBASE_HPP

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    template<typename KernelImplementation>
    class DMSystemMatrixVectorizedIdentityMPIBase : public SGPP::datadriven::DMSystemMatrixBase {
      protected:
        VectorizationType vecMode_;
        size_t numTrainingInstances_;
        size_t numPatchedTrainingInstances_;

        /// Pointer to the grid's gridstorage object
        SGPP::base::GridStorage* storage_;
        /// Member to store the sparse grid's levels for better vectorization
        SGPP::base::DataMatrix* level_;
        /// Member to store the sparse grid's indices for better vectorization
        SGPP::base::DataMatrix* index_;
        /// Member to store for masks per grid point for better vectorization of modlinear operations
        SGPP::base::DataMatrix* mask_;
        /// Member to store offsets per grid point for better vecotrization of modlinear operations
        SGPP::base::DataMatrix* offset_;

        KernelImplementation kernel_;

        /// only allocate temporary arrays once
        SGPP::base::DataVector* tempData;
        SGPP::base::DataVector* result_tmp;
      public:
        DMSystemMatrixVectorizedIdentityMPIBase(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
          : DMSystemMatrixBase(trainData, lambda),
            vecMode_(vecMode),
            numTrainingInstances_(0),
            numPatchedTrainingInstances_(0),
            storage_(SparseGrid.getStorage()),
            level_(NULL),
            index_(NULL),
            mask_(NULL),
            offset_(NULL) {
          // handle unsupported vector extensions
          if (this->vecMode_ != X86SIMD &&
              this->vecMode_ != MIC &&
              this->vecMode_ != Hybrid_X86SIMD_MIC &&
              this->vecMode_ != OpenCL &&
              this->vecMode_ != ArBB &&
              this->vecMode_ != Hybrid_X86SIMD_OpenCL) {
            throw SGPP::base::operation_exception("DMSystemMatrixVectorizedIdentityAllreduce : un-supported vector extension!");
          }

          this->dataset_ = new SGPP::base::DataMatrix(trainData);
          this->numTrainingInstances_ = this->dataset_->getNrows();
          this->numPatchedTrainingInstances_ = SGPP::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);
          std::cout << "Padding Dataset to " << numPatchedTrainingInstances_ << " Instances. " << std::endl;
          this->tempData = new SGPP::base::DataVector(this->numPatchedTrainingInstances_);

          if (this->vecMode_ != ArBB) {
            this->dataset_->transpose();
          }

          this->result_tmp = new SGPP::base::DataVector(storage_->size());
        }

        virtual ~DMSystemMatrixVectorizedIdentityMPIBase() {
          if (this->level_ != NULL)
            delete this->level_;

          if (this->index_ != NULL)
            delete this->index_;

          if (this->mask_ != NULL)
            delete this->mask_;

          if (this->offset_ != NULL)
            delete this->offset_;

          delete this->result_tmp;
          delete this->tempData;
          delete this->dataset_;
        }

        virtual void rebuildLevelAndIndex() {
          LevelIndexMaskOffsetHelper::rebuild<KernelImplementation::kernelType, DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation> >(this);

          if (this->result_tmp != NULL) {
            delete this->result_tmp;
          }

          this->result_tmp = new SGPP::base::DataVector(storage_->size());
          this->kernel_.resetKernel();
        }
        friend struct LevelIndexMaskOffsetHelper::rebuild<KernelImplementation::kernelType, DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation> >;
    };

  }
}

#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYMPIBASE_HPP