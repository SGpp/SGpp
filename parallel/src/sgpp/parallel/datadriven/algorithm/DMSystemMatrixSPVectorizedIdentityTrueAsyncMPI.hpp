// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITYTRUEASYNCMPI_H
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITYTRUEASYNCMPI_H

#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPIBase.hpp>
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/globaldef.hpp>

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace parallel {

template <typename KernelImplementation>
class DMSystemMatrixSPVectorizedIdentityTrueAsyncMPI
    : public sgpp::parallel::DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation> {
 private:
  /// how to distribute storage array across processes
  int* _mpi_grid_sizes;
  int* _mpi_grid_offsets;

  /// how to distribute dataset across processes
  int* _mpi_data_sizes;
  int* _mpi_data_offsets;

 public:
  /**
   * Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to sgpp::base::DataMatrix that contains the training data
   * @param lambda the lambda, the regression parameter
   * @param vecMode vectorization mode
   */
  DMSystemMatrixSPVectorizedIdentityTrueAsyncMPI(sgpp::base::Grid& SparseGrid,
                                                 sgpp::base::DataMatrixSP& trainData, float lambda,
                                                 VectorizationType vecMode)
      : DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation>(SparseGrid, trainData,
                                                                        lambda, vecMode) {
    size_t mpi_size = sgpp::parallel::myGlobalMPIComm->getNumRanks();

    // arrays for distribution settings
    _mpi_data_sizes = new int[mpi_size];
    _mpi_data_offsets = new int[mpi_size];

    sgpp::parallel::PartitioningTool::calcDistribution(
        this->numPatchedTrainingInstances_, mpi_size, _mpi_data_sizes, _mpi_data_offsets,
        sgpp::parallel::DMVectorizationPaddingAssistant::getVecWidthSP(this->vecMode_));

    // arrays for distribution settings
    _mpi_grid_sizes = new int[mpi_size];
    _mpi_grid_offsets = new int[mpi_size];

    rebuildLevelAndIndex();

    // mult: distribute calculations over dataset
    // multTranspose: distribute calculations over grid
  }

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixSPVectorizedIdentityTrueAsyncMPI() {
    delete[] this->_mpi_grid_sizes;
    delete[] this->_mpi_grid_offsets;
    delete[] this->_mpi_data_sizes;
    delete[] this->_mpi_data_offsets;
  }

  virtual void mult(sgpp::base::DataVectorSP& alpha, sgpp::base::DataVectorSP& result) {
#ifdef X86_MIC_SYMMETRIC
    myGlobalMPIComm->broadcastSPGridCoefficientsFromRank0(alpha);
#endif
    result.setAll(0.0f);
    this->tempData->setAll(0.0f);
    float* ptrTemp = this->tempData->getPointer();

    size_t mpi_size = sgpp::parallel::myGlobalMPIComm->getNumRanks();
    size_t mpi_myrank = sgpp::parallel::myGlobalMPIComm->getMyRank();

    MPI_Request* dataSendReqs = new MPI_Request[mpi_size];
    MPI_Request* dataRecvReqs = new MPI_Request[mpi_size];  // allocating a little more than
                                                            // necessary, otherwise complicated
                                                            // index computations needed
    int* tagsData = new int[mpi_size];

    for (size_t i = 0; i < mpi_size; i++) {
      tagsData[i] = _mpi_data_offsets[i] * 2 + 2;
    }

    myGlobalMPIComm->IrecvFromAllSP(ptrTemp, 1, _mpi_data_sizes, _mpi_data_offsets, tagsData,
                                    dataRecvReqs);

    this->myTimer_->start();

    size_t dataProcessChunkStart = _mpi_data_offsets[mpi_myrank];
    size_t dataProcessChunkEnd = dataProcessChunkStart + _mpi_data_sizes[mpi_myrank];
    size_t gridProcessChunkStart = _mpi_grid_offsets[mpi_myrank];
    size_t gridProcessChunkEnd = gridProcessChunkStart + _mpi_grid_sizes[mpi_myrank];
    int idx;
#pragma omp parallel
    {
      this->kernel_.mult(this->level_, this->index_, this->mask_, this->offset_, this->dataset_,
                         alpha, *(this->tempData), 0, alpha.getSize(), dataProcessChunkStart,
                         dataProcessChunkEnd);

#pragma omp barrier
// make sure that all threads finished their part, so that we can
// safely overwrite the padded range and send the results to all other procs
#pragma omp single
      {
        // patch result -> set additional entries zero
        for (size_t i = this->numTrainingInstances_; i < this->numPatchedTrainingInstances_; i++) {
          this->tempData->set(i, 0.0f);
        }
      }
// implicit openMP barrier here (after omp single): we have to make sure that we do not
// read non-zeros from padded temp range

#pragma omp master
      {
        sgpp::parallel::myGlobalMPIComm->IsendToAllSP(&ptrTemp[dataProcessChunkStart],
                                                      _mpi_data_sizes[mpi_myrank],
                                                      tagsData[mpi_myrank], dataSendReqs);
      }
      // we don't need to wait for the sendreqs to finish because we only read from temp after its
      // calculation

      // while the data is transfered we already calculate multTrans with this part of the grid and
      // the temp-values computed by this process
      this->kernel_.multTranspose(this->level_, this->index_, this->mask_, this->offset_,
                                  this->dataset_, *(this->tempData), result, gridProcessChunkStart,
                                  gridProcessChunkEnd, dataProcessChunkStart, dataProcessChunkEnd);

      // after this, we receive the temp chunks from all the other processes and do the calculations
      // for them
      while (true) {
#pragma omp master
        {
          myGlobalMPIComm->waitForAnyRequest(sgpp::parallel::myGlobalMPIComm->getNumRanks(),
                                             dataRecvReqs, &idx);
        }
#pragma omp barrier
        // let all threads wait here

        if (idx == MPI_UNDEFINED) {
          // no more active request, everything is done
          break;
        }

        size_t dataChunkStart = _mpi_data_offsets[idx];
        size_t dataChunkEnd = dataChunkStart + _mpi_data_sizes[idx];

// make sure that all threads have read the same idx
// if one thread finished computation and restarts the loop and reads the next index before another
// thread could read the idx, the computation would be erroneous (two times computation of
// same index or not computing last index)
#pragma omp barrier

        this->kernel_.multTranspose(this->level_, this->index_, this->mask_, this->offset_,
                                    this->dataset_, *(this->tempData), result,
                                    gridProcessChunkStart, gridProcessChunkEnd, dataChunkStart,
                                    dataChunkEnd);
      }
    }
    // send result of this process to all other processes

    this->computeTimeMultTrans_ += this->myTimer_->stop();
    sgpp::parallel::myGlobalMPIComm->dataVectorSPAllToAll(result, _mpi_grid_offsets,
                                                          _mpi_grid_sizes);
    this->completeTimeMultTrans_ += this->myTimer_->stop();

    result.axpy(static_cast<float>(this->numTrainingInstances_) * this->lambda_, alpha);

    if (mpi_myrank == 0) std::cout << "*";

    delete[] dataSendReqs;
    delete[] dataRecvReqs;
    delete[] tagsData;
  }

  virtual void generateb(sgpp::base::DataVectorSP& classes, sgpp::base::DataVectorSP& b) {
    size_t mpi_myrank = sgpp::parallel::myGlobalMPIComm->getMyRank();
    b.setAll(0.0f);

    this->tempData->setAll(0.0);
    this->tempData->copyFrom(classes);

    this->myTimer_->start();
#pragma omp parallel
    {
      this->kernel_.multTranspose(this->level_, this->index_, this->mask_, this->offset_,
                                  this->dataset_, *(this->tempData), b,
                                  _mpi_grid_offsets[mpi_myrank],
                                  _mpi_grid_offsets[mpi_myrank] + _mpi_grid_sizes[mpi_myrank], 0,
                                  this->numPatchedTrainingInstances_);
    }
    this->computeTimeMultTrans_ += this->myTimer_->stop();
    sgpp::parallel::myGlobalMPIComm->dataVectorSPAllToAll(b, _mpi_grid_offsets, _mpi_grid_sizes);
    this->completeTimeMultTrans_ += this->myTimer_->stop();
  }

  virtual void rebuildLevelAndIndex() {
    DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation>::rebuildLevelAndIndex();

    size_t mpi_size = sgpp::parallel::myGlobalMPIComm->getNumRanks();

    sgpp::parallel::PartitioningTool::calcDistribution(this->storage_.getSize(), mpi_size,
                                                       _mpi_grid_sizes, _mpi_grid_offsets, 1);
  }
};

}  // namespace parallel
}  // namespace sgpp

#endif  // DMSYSTEMMATRIXSPVECTORIZEDIDENTITYTRUEASYNCMPI_H
