/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITYMPIBASE_HPP
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITYMPIBASE_HPP

#include <sgpp/datadriven/algorithm/DMSystemMatrixBaseSP.hpp>
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp>


namespace sg {
  namespace parallel {
    template<typename KernelImplementation>
    class DMSystemMatrixSPVectorizedIdentityMPIBase : public sg::datadriven::DMSystemMatrixBaseSP {
      protected:
        VectorizationType vecMode_;
        size_t numTrainingInstances_;
        size_t numPatchedTrainingInstances_;

        /// Pointer to the grid's gridstorage object
        sg::base::GridStorage* storage_;
        /// Member to store the sparse grid's levels for better vectorization
        sg::base::DataMatrixSP* level_;
        /// Member to store the sparse grid's indices for better vectorization
        sg::base::DataMatrixSP* index_;
        /// Member to store for masks per grid point for better vectorization of modlinear operations
        sg::base::DataMatrixSP* mask_;
        /// Member to store offsets per grid point for better vecotrization of modlinear operations
        sg::base::DataMatrixSP* offset_;

        KernelImplementation kernel_;

        /// only allocate temporary arrays once
        sg::base::DataVectorSP* tempData;
        sg::base::DataVectorSP* result_tmp;
      public:
        DMSystemMatrixSPVectorizedIdentityMPIBase(sg::base::Grid& SparseGrid, sg::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode)
          : DMSystemMatrixBaseSP(trainData, lambda),
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
            throw sg::base::operation_exception("DMSystemMatrixVectorizedIdentityAllreduce : un-supported vector extension!");
          }

          this->dataset_ = new sg::base::DataMatrixSP(trainData);
          this->numTrainingInstances_ = this->dataset_->getNrows();
          this->numPatchedTrainingInstances_ = sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);
          std::cout << "Padding Dataset to " << numPatchedTrainingInstances_ << " Instances. " << std::endl;
          this->tempData = new sg::base::DataVectorSP(this->numPatchedTrainingInstances_);

          if (this->vecMode_ != ArBB) {
            this->dataset_->transpose();
          }

          this->result_tmp = new sg::base::DataVectorSP(storage_->size());
        }

        virtual ~DMSystemMatrixSPVectorizedIdentityMPIBase() {
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
          LevelIndexMaskOffsetHelperSP::rebuild<KernelImplementation::kernelType, DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation> >(this);

          if (this->result_tmp != NULL) {
            delete this->result_tmp;
          }

          this->result_tmp = new sg::base::DataVectorSP(storage_->size());
          kernel_.resetKernel();
        }
        friend struct LevelIndexMaskOffsetHelperSP::rebuild<KernelImplementation::kernelType, DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation> >;
    };

  }
}

#endif // DMSYSTEMMATRIXSPVECTORIZEDIDENTITYMPIBASE_HPP
