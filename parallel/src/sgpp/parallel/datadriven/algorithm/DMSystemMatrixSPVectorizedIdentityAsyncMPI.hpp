// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITYASYNCMPI_HPP
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITYASYNCMPI_HPP

#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPIBase.hpp>
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <strstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/globaldef.hpp>


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
    template<typename KernelImplementation>
    class DMSystemMatrixSPVectorizedIdentityAsyncMPI : public SGPP::parallel::DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation> {
      public:
        /**
         * Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to SGPP::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         * @param vecMode vectorization mode
         */
        DMSystemMatrixSPVectorizedIdentityAsyncMPI(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode)
          : DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation>(SparseGrid, trainData, lambda, vecMode) {
          size_t mpi_size = SGPP::parallel::myGlobalMPIComm->getNumRanks();

          /* calculate distribution of data */
          _chunkCountPerProcData = 2;
          _mpi_data_sizes = new int[_chunkCountPerProcData * mpi_size];
          _mpi_data_offsets = new int[_chunkCountPerProcData * mpi_size];
          PartitioningTool::calcMPIChunkedDistribution(this->numPatchedTrainingInstances_, _chunkCountPerProcData, _mpi_data_sizes, _mpi_data_offsets, SGPP::parallel::DMVectorizationPaddingAssistant::getVecWidthSP(this->vecMode_));

          if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
            std::cout << "Max size per chunk Data: " << _mpi_data_sizes[0] << std::endl;
          }

          _mpi_grid_sizes = NULL; // allocation in rebuildLevelAndIndex();
          _mpi_grid_offsets = NULL; // allocation in rebuildLevelAndIndex();
          rebuildLevelAndIndex();

          // mult: distribute calculations over dataset
          // multTranspose: distribute calculations over grid
        }

        /**
         * Destructor
         */
        virtual ~DMSystemMatrixSPVectorizedIdentityAsyncMPI() {
          delete[] this->_mpi_grid_sizes;
          delete[] this->_mpi_grid_offsets;
          delete[] this->_mpi_data_sizes;
          delete[] this->_mpi_data_offsets;
        }

        virtual void mult(SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& result) {
#ifdef X86_MIC_SYMMETRIC
          myGlobalMPIComm->broadcastSPGridCoefficientsFromRank0(alpha);
#endif
          SGPP::base::DataVectorSP temp(this->numPatchedTrainingInstances_);
          result.setAll(0.0f);
          temp.setAll(0.0f);
          float* ptrResult = result.getPointer();
          float* ptrTemp = temp.getPointer();

          size_t mpi_size = SGPP::parallel::myGlobalMPIComm->getNumRanks();
          size_t mpi_myrank = SGPP::parallel::myGlobalMPIComm->getMyRank();

          size_t totalChunkCountGrid = _chunkCountPerProcGrid * mpi_size;
          size_t totalChunkCountData = _chunkCountPerProcData * mpi_size;

          /* setup MPI_Requests, tags and post receives for data */
          MPI_Request* dataRecvReqs = new MPI_Request[totalChunkCountData]; //allocating a little more than necessary, otherwise complicated index computations needed
          int* tagsData = new int[totalChunkCountData];

          for (size_t i = 0; i < totalChunkCountData; i++) {
            tagsData[i] = (int)(i * 2 + 2);
          }

          SGPP::parallel::myGlobalMPIComm->IrecvFromAllSP(ptrTemp, _chunkCountPerProcData, _mpi_data_sizes, _mpi_data_offsets, tagsData, dataRecvReqs);

          /* setup MPI_Requests, tags and post receives for grid */
          MPI_Request* gridRecvReqs = new MPI_Request[totalChunkCountGrid]; //allocating a little more than necessary, otherwise complicated index computations needed
          int* tagsGrid = new int[totalChunkCountGrid];

          for (size_t i = 0; i < totalChunkCountGrid; i++) {
            tagsGrid[i] = (int)(i * 2 + 3);
          }

          SGPP::parallel::myGlobalMPIComm->IrecvFromAllSP(ptrResult, _chunkCountPerProcGrid, _mpi_grid_sizes, _mpi_grid_offsets, tagsGrid, gridRecvReqs);
          MPI_Request* dataSendReqs = new MPI_Request[totalChunkCountData];
          MPI_Request* gridSendReqs = new MPI_Request[totalChunkCountGrid];

          this->myTimer_->start();
          #pragma omp parallel
          {
            size_t myDataChunkStart = mpi_myrank * _chunkCountPerProcData;
            size_t myDataChunkEnd = (mpi_myrank + 1) * _chunkCountPerProcData;

            for (size_t chunkIndex = myDataChunkStart; chunkIndex < myDataChunkEnd; chunkIndex++) {
              size_t start = _mpi_data_offsets[chunkIndex];
              size_t end =  start + _mpi_data_sizes[chunkIndex];
              this->kernel_.mult(
                this->level_,
                this->index_,
                this->mask_,
                this->offset_,
                this->dataset_,
                alpha,
                temp,
                0,
                alpha.getSize(),
                start,
                end);
              #pragma omp barrier
              #pragma omp master // the non-sending processes can already continue with execution
              {
                myGlobalMPIComm->IsendToAllSP(&ptrTemp[start], _mpi_data_sizes[chunkIndex], tagsData[chunkIndex], &dataSendReqs[(chunkIndex - myDataChunkStart)*mpi_size]);
              }
            }

            #pragma omp single
            {
              // patch result -> set additional entries zero
              // only done for processes that need this part of the temp data for multTrans
              for (size_t i = this->numTrainingInstances_; i < this->numPatchedTrainingInstances_; i++) {
                temp.set(i, 0.0f);
              }
            }

            #pragma omp master
            {
              double computationTime = this->myTimer_->stop();
              this->computeTimeMult_ += computationTime;
              myGlobalMPIComm->waitForAllRequests(totalChunkCountData, dataRecvReqs);

              // we don't really need to wait for the sends to
              // finish as we don't need (in particular not modify) temp
              // advantage: it's faster like this
              // myGlobalMPIComm->waitForAllRequests(totalChunkCountData, dataSendReqs);
              this->completeTimeMult_ += this->myTimer_->stop();

              this->myTimer_->start();

            }
            #pragma omp barrier

            size_t myGridChunkStart = mpi_myrank * _chunkCountPerProcGrid;
            size_t myGridChunkEnd = (mpi_myrank + 1) * _chunkCountPerProcGrid;

            for (size_t chunkIndex = myGridChunkStart; chunkIndex < myGridChunkEnd; chunkIndex++) {
              size_t start = _mpi_grid_offsets[chunkIndex];
              size_t end =  start + _mpi_grid_sizes[chunkIndex];
              this->kernel_.multTranspose(
                this->level_,
                this->index_,
                this->mask_,
                this->offset_,
                this->dataset_,
                temp,
                result,
                start,
                end,
                0,
                this->numPatchedTrainingInstances_
              );
              #pragma omp barrier
              #pragma omp master // the non-sending processes can already continue with execution
              {
                myGlobalMPIComm->IsendToAllSP(&ptrResult[start], _mpi_grid_sizes[chunkIndex], tagsGrid[chunkIndex], &gridSendReqs[(chunkIndex - myGridChunkStart)*mpi_size]);
              }
            }
          }
          this->computeTimeMultTrans_ += this->myTimer_->stop();
          myGlobalMPIComm->waitForAllRequests(totalChunkCountGrid, gridRecvReqs);
          myGlobalMPIComm->waitForAllRequests(totalChunkCountGrid, gridSendReqs);

          this->completeTimeMultTrans_ += this->myTimer_->stop();

          result.axpy(static_cast<float>(this->numTrainingInstances_)*this->lambda_, alpha);

          if (mpi_myrank == 0) std::cout << "*";

          delete[] dataSendReqs;
          delete[] gridSendReqs;
          delete[] dataRecvReqs;
          delete[] gridRecvReqs;
          delete[] tagsData;
          delete[] tagsGrid;
        } //end mult

        virtual void generateb(SGPP::base::DataVectorSP& classes, SGPP::base::DataVectorSP& b) {
          size_t mpi_size = SGPP::parallel::myGlobalMPIComm->getNumRanks();
          size_t mpi_myrank = SGPP::parallel::myGlobalMPIComm->getMyRank();

          float* ptrB = b.getPointer();
          b.setAll(0.0f);

          SGPP::base::DataVectorSP myClasses(classes);

          // Apply padding
          if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
            myClasses.resizeZero(this->numPatchedTrainingInstances_);
          }

          size_t totalChunkCount = mpi_size * _chunkCountPerProcGrid;
          MPI_Request* gridRecvReqs = new MPI_Request[totalChunkCount]; //allocating a little more than necessary, otherwise complicated index computations needed
          int* tags = new int[totalChunkCount];

          for (size_t i = 0; i < totalChunkCount; i++) {
            tags[i] = (int)(i + 1);
          }

          SGPP::parallel::myGlobalMPIComm->IrecvFromAllSP(ptrB, _chunkCountPerProcGrid, _mpi_grid_sizes, _mpi_grid_offsets, tags, gridRecvReqs);
          MPI_Request* gridSendReqs = new MPI_Request[totalChunkCount];
          this->myTimer_->start();
          #pragma omp parallel
          {
            size_t myGridChunkStart = mpi_myrank * _chunkCountPerProcGrid;
            size_t myGridChunkEnd = (mpi_myrank + 1) * _chunkCountPerProcGrid;

            for (size_t chunkIndex = myGridChunkStart; chunkIndex < myGridChunkEnd; chunkIndex++) {
              size_t start = _mpi_grid_offsets[chunkIndex];
              size_t end =  start + _mpi_grid_sizes[chunkIndex];
              this->kernel_.multTranspose(
                this->level_,
                this->index_,
                this->mask_,
                this->offset_,
                this->dataset_,
                myClasses,
                b,
                start,
                end,
                0,
                this->numPatchedTrainingInstances_
              );
              #pragma omp barrier
              #pragma omp master // the non-sending processes can already continue with execution
              {
                myGlobalMPIComm->IsendToAllSP(&ptrB[start], _mpi_grid_sizes[chunkIndex], tags[chunkIndex], &gridSendReqs[(chunkIndex - myGridChunkStart)*mpi_size]);
              }
            }
          }
          this->computeTimeMultTrans_ += this->myTimer_->stop();
          myGlobalMPIComm->waitForAllRequests(totalChunkCount, gridRecvReqs);
          myGlobalMPIComm->waitForAllRequests(totalChunkCount, gridSendReqs);
          this->completeTimeMultTrans_ += this->myTimer_->stop();
          delete[] gridRecvReqs;
          delete[] gridSendReqs;
          delete[] tags;
        }

        virtual void rebuildLevelAndIndex() {
          DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation>::rebuildLevelAndIndex();

          if (_mpi_grid_sizes != NULL) {
            delete[] _mpi_grid_sizes;
          }

          if (_mpi_grid_offsets != NULL) {
            delete[] _mpi_grid_offsets;
          }

          size_t mpi_size = myGlobalMPIComm->getNumRanks();
          size_t sendChunkSize = 16;
          size_t sizePerProc = this->storage_->size() / mpi_size;
          std::max<size_t>(sizePerProc / sendChunkSize, 1);

          _chunkCountPerProcGrid = 1;

          if (myGlobalMPIComm->getMyRank() == 0) {
            std::cout << "chunksperproc grid: " << _chunkCountPerProcGrid << "; total # chunks: " << _chunkCountPerProcGrid* mpi_size << std::endl;
          }

          _mpi_grid_sizes = new int[_chunkCountPerProcGrid * mpi_size];
          _mpi_grid_offsets = new int[_chunkCountPerProcGrid * mpi_size];
          PartitioningTool::calcMPIChunkedDistribution(this->storage_->size(), _chunkCountPerProcGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
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
    };

  }
}
#endif /* DMSYSTEMMATRIXSPVECTORIZEDIDENTITYASYNCMPI_HPP */