// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP

#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPIBase.hpp>
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/globaldef.hpp>

#include <strstream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SGPP {
namespace parallel {

/**
 * Class that implements the virtual class SGPP::base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
template <typename KernelImplementation>
class DMSystemMatrixVectorizedIdentityOnesidedMPI
    : public SGPP::parallel::DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation> {
 public:
  /**
   * Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to SGPP::base::DataMatrix that contains the training data
   * @param lambda the lambda, the regression parameter
   * @param vecMode vectorization mode
   */
  DMSystemMatrixVectorizedIdentityOnesidedMPI(SGPP::base::Grid& SparseGrid,
                                              SGPP::base::DataMatrix& trainData, double lambda,
                                              VectorizationType vecMode)
      : DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation>(SparseGrid, trainData, lambda,
                                                                      vecMode),
        _mpi_grid_window_buffer(NULL),
        _mpi_data_window_buffer(NULL),
        _mpi_grid_window(),
        _mpi_data_window() {
    size_t mpi_size = SGPP::parallel::myGlobalMPIComm->getNumRanks();

    /* calculate distribution of data */
    _chunkCountPerProcData = 2;
    _mpi_data_sizes = new int[_chunkCountPerProcData * mpi_size];
    _mpi_data_offsets = new int[_chunkCountPerProcData * mpi_size];
    PartitioningTool::calcMPIChunkedDistribution(
        this->numPatchedTrainingInstances_, _chunkCountPerProcData, _mpi_data_sizes,
        _mpi_data_offsets,
        SGPP::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

    if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
      std::cout << "Max size per chunk Data: " << _mpi_data_sizes[0] << std::endl;
    }

    _mpi_grid_sizes = NULL;    // allocation in rebuildLevelAndIndex();
    _mpi_grid_offsets = NULL;  // allocation in rebuildLevelAndIndex();
    rebuildLevelAndIndex();

    createAndInitDataBuffers();
    // mult: distribute calculations over dataset
    // multTranspose: distribute calculations over grid
  }

  void createAndInitGridBuffers() {
    if (_mpi_grid_window_buffer != NULL) {
      delete _mpi_grid_window_buffer;
      MPI_Win_free(&_mpi_grid_window);
    }

    _mpi_grid_window_buffer = new SGPP::base::DataVector(this->storage_.getSize());

    MPI_Win_create(_mpi_grid_window_buffer->getPointer(), this->storage_.getSize() * sizeof(double),
                   static_cast<int>(sizeof(double)), MPI_INFO_NULL, MPI_COMM_WORLD,
                   &_mpi_grid_window);
    MPI_Win_fence(MPI_MODE_NOSUCCEED, _mpi_grid_window);
  }

  void createAndInitDataBuffers() {
    // create MPI windows for databuffer
    if (_mpi_data_window_buffer != NULL) {
      delete _mpi_data_window_buffer;
      MPI_Win_free(&_mpi_data_window);
    }

    _mpi_data_window_buffer = new SGPP::base::DataVector(this->numPatchedTrainingInstances_);
    MPI_Win_create(
        _mpi_data_window_buffer->getPointer(), this->numPatchedTrainingInstances_ * sizeof(double),
        static_cast<int>(sizeof(double)), MPI_INFO_NULL, MPI_COMM_WORLD, &_mpi_data_window);
    MPI_Win_fence(MPI_MODE_NOSUCCEED, _mpi_data_window);
  }

  /**
   * Destructor
   */
  virtual ~DMSystemMatrixVectorizedIdentityOnesidedMPI() {
    delete[] this->_mpi_grid_sizes;
    delete[] this->_mpi_grid_offsets;
    delete[] this->_mpi_data_sizes;
    delete[] this->_mpi_data_offsets;
    MPI_Win_free(&_mpi_data_window);
    MPI_Win_free(&_mpi_grid_window);
  }

  virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
#ifdef X86_MIC_SYMMETRIC
    myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);
#endif
    SGPP::base::DataVector temp(this->numPatchedTrainingInstances_);

    result.setAll(0.0);
    temp.setAll(0.0);
    _mpi_grid_window_buffer->setAll(0.0);
    _mpi_data_window_buffer->setAll(0.0);
    double* ptrResult = result.getPointer();
    double* ptrTemp = temp.getPointer();

    size_t mpi_myrank = SGPP::parallel::myGlobalMPIComm->getMyRank();

    // we expect that there were no RMA calls before this line
    MPI_Win_fence(MPI_MODE_NOPRECEDE, _mpi_grid_window);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, _mpi_data_window);

    this->myTimer_->start();
#pragma omp parallel
    {
      size_t myDataChunkStart = mpi_myrank * _chunkCountPerProcData;
      size_t myDataChunkEnd = (mpi_myrank + 1) * _chunkCountPerProcData;

      for (size_t chunkIndex = myDataChunkStart; chunkIndex < myDataChunkEnd; chunkIndex++) {
        size_t start = _mpi_data_offsets[chunkIndex];
        size_t end = start + _mpi_data_sizes[chunkIndex];
        this->kernel_.mult(this->level_, this->index_, this->mask_, this->offset_, this->dataset_,
                           alpha, temp, 0, alpha.getSize(), start, end);
#pragma omp barrier
#pragma omp master  // the non-sending processes can already continue with execution
        {
          myGlobalMPIComm->putToAll(&ptrTemp[start], start, _mpi_data_sizes[chunkIndex],
                                    _mpi_data_window);
        }
      }

// this barrier is absolutely necessary, as it guarantees that the last
// put-call has been made. This is required as a MPI_Win_fence
// follows. If that fence would be called before the master thread has finished,
// the last put won't be taken into account for the synchronization.
#pragma omp barrier

#pragma omp master
      {
        this->computeTimeMult_ += this->myTimer_->stop();

        // we expect that there is not put to the buffer after this line
        // and we did no store to the buffer locally (only via put to the same proc)
        // and there won't be RMA calls until the next fence
        MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED, _mpi_data_window);

        // patch result -> set additional entries zero
        for (size_t i = this->numTrainingInstances_; i < this->numPatchedTrainingInstances_; i++) {
          _mpi_data_window_buffer->set(i, 0.0f);
        }

        this->completeTimeMult_ += this->myTimer_->stop();
        this->myTimer_->start();
      }
#pragma omp barrier

      size_t myGridChunkStart = mpi_myrank * _chunkCountPerProcGrid;
      size_t myGridChunkEnd = (mpi_myrank + 1) * _chunkCountPerProcGrid;

      for (size_t chunkIndex = myGridChunkStart; chunkIndex < myGridChunkEnd; chunkIndex++) {
        size_t start = _mpi_grid_offsets[chunkIndex];
        size_t end = start + _mpi_grid_sizes[chunkIndex];
        this->kernel_.multTranspose(this->level_, this->index_, this->mask_, this->offset_,
                                    this->dataset_, *_mpi_data_window_buffer, result, start, end, 0,
                                    this->numPatchedTrainingInstances_);
#pragma omp barrier
#pragma omp master  // the non-sending processes can already continue with execution
        {
          myGlobalMPIComm->putToAll(&ptrResult[start], start, _mpi_grid_sizes[chunkIndex],
                                    _mpi_grid_window);
        }
      }
    }
    this->computeTimeMultTrans_ += this->myTimer_->stop();

    MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED, _mpi_grid_window);
    this->completeTimeMultTrans_ += this->myTimer_->stop();
    result.copyFrom(*_mpi_grid_window_buffer);
    result.axpy(static_cast<double>(this->numTrainingInstances_) * this->lambda_, alpha);
  }  // end mult

  virtual void generateb(SGPP::base::DataVector& classes, SGPP::base::DataVector& b) {
    size_t mpi_myrank = SGPP::parallel::myGlobalMPIComm->getMyRank();

    _mpi_grid_window_buffer->setAll(0.0);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, _mpi_grid_window);

    double* ptrB = b.getPointer();
    b.setAll(0.0);

    SGPP::base::DataVector myClasses(classes);

    // Apply padding
    if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
      myClasses.resizeZero(this->numPatchedTrainingInstances_);
    }

    this->myTimer_->start();
#pragma omp parallel
    {
      size_t myGridChunkStart = mpi_myrank * _chunkCountPerProcGrid;
      size_t myGridChunkEnd = (mpi_myrank + 1) * _chunkCountPerProcGrid;

      for (size_t chunkIndex = myGridChunkStart; chunkIndex < myGridChunkEnd; chunkIndex++) {
        size_t start = _mpi_grid_offsets[chunkIndex];
        size_t end = start + _mpi_grid_sizes[chunkIndex];
        this->kernel_.multTranspose(this->level_, this->index_, this->mask_, this->offset_,
                                    this->dataset_, myClasses, b, start, end, 0,
                                    this->numPatchedTrainingInstances_);
#pragma omp barrier
#pragma omp master  // the non-sending processes can already continue with execution
        {
          myGlobalMPIComm->putToAll(&ptrB[start], start, _mpi_grid_sizes[chunkIndex],
                                    _mpi_grid_window);
        }
      }
    }

    this->computeTimeMultTrans_ += this->myTimer_->stop();
    MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED, _mpi_grid_window);
    this->completeTimeMultTrans_ += this->myTimer_->stop();
    b.copyFrom(*_mpi_grid_window_buffer);
  }

  virtual void rebuildLevelAndIndex() {
    DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation>::rebuildLevelAndIndex();

    if (_mpi_grid_sizes != NULL) {
      delete[] _mpi_grid_sizes;
    }

    if (_mpi_grid_offsets != NULL) {
      delete[] _mpi_grid_offsets;
    }

    size_t mpi_size = myGlobalMPIComm->getNumRanks();

    size_t sendChunkSize = 2;
    size_t sizePerProc = this->storage_.getSize() / mpi_size;
    std::max<size_t>(sizePerProc / sendChunkSize, 1);

    _chunkCountPerProcGrid = 1;

    if (myGlobalMPIComm->getMyRank() == 0) {
      std::cout << "chunksperproc grid: " << _chunkCountPerProcGrid
                << "; total # chunks: " << _chunkCountPerProcGrid * mpi_size << std::endl;
    }

    _mpi_grid_sizes = new int[_chunkCountPerProcGrid * mpi_size];
    _mpi_grid_offsets = new int[_chunkCountPerProcGrid * mpi_size];
    PartitioningTool::calcMPIChunkedDistribution(this->storage_.getSize(), _chunkCountPerProcGrid,
                                                 _mpi_grid_sizes, _mpi_grid_offsets, 1);
    createAndInitGridBuffers();
  }

 private:
  /// how to distribute storage array across processes
  int* _mpi_grid_sizes;
  int* _mpi_grid_offsets;

  /// how to distribute dataset across processes
  int* _mpi_data_sizes;
  int* _mpi_data_offsets;

  /// into how many chunks should data and grid be partitioned
  size_t _chunkCountPerProcData;
  size_t _chunkCountPerProcGrid;

  /// MPI windows
  SGPP::base::DataVector* _mpi_grid_window_buffer;
  SGPP::base::DataVector* _mpi_data_window_buffer;
  MPI_Win _mpi_grid_window;
  MPI_Win _mpi_data_window;
};
}  // namespace parallel
}  // namespace SGPP
#endif  // DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP
