/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPI.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg {
  namespace parallel {

    DMSystemMatrixVectorizedIdentityMPI::DMSystemMatrixVectorizedIdentityMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
      : DMSystemMatrixBase(trainData, lambda), vecMode_(vecMode), numTrainingInstances_(0), numPatchedTrainingInstances_(0), m_grid(SparseGrid) {
      // handle unsupported vector extensions
      if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC && this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL) {
        throw sg::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
      }

      this->dataset_ = new sg::base::DataMatrix(trainData);
      this->numTrainingInstances_ = this->dataset_->getNrows();
      this->numPatchedTrainingInstances_ = sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);

      if (this->vecMode_ != ArBB) {
        this->dataset_->transpose();
      }

      size_t mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
      size_t mpi_rank = sg::parallel::myGlobalMPIComm->getMyRank();

      // arrays for distribution settings
      _mpi_grid_sizes = new int[mpi_size];
      _mpi_grid_offsets = new int[mpi_size];
      _mpi_data_sizes = new int[mpi_size];
      _mpi_data_offsets = new int[mpi_size];

      // calculate distribution
      sg::parallel::PartitioningTool::calcDistribution(m_grid.getStorage()->size(), mpi_size, _mpi_grid_sizes, _mpi_grid_offsets, 1);
      sg::parallel::PartitioningTool::calcDistribution(this->numPatchedTrainingInstances_, mpi_size, _mpi_data_sizes, _mpi_data_offsets,
          sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

      this->B_ = sg::op_factory::createOperationMultipleEvalVectorized(m_grid, this->vecMode_, this->dataset_,
                 _mpi_grid_offsets[mpi_rank],
                 _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank],
                 _mpi_data_offsets[mpi_rank],
                 _mpi_data_offsets[mpi_rank] + _mpi_data_sizes[mpi_rank]
                                                                      );
      waitting_time = 0.0;
    }

    DMSystemMatrixVectorizedIdentityMPI::~DMSystemMatrixVectorizedIdentityMPI() {
      delete this->B_;
      delete this->dataset_;

      delete[] this->_mpi_grid_sizes;
      delete[] this->_mpi_grid_offsets;
      delete[] this->_mpi_data_sizes;
      delete[] this->_mpi_data_offsets;
      double* allwaittime = new double [sg::parallel::myGlobalMPIComm->getNumRanks()];
      double* allcomptimemultTrans = new double [sg::parallel::myGlobalMPIComm->getNumRanks()];
      double* allcomptimemult = new double [sg::parallel::myGlobalMPIComm->getNumRanks()];
      MPI_Gather(&waitting_time, 1, MPI_DOUBLE, allwaittime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&computeTimeMultTrans_, 1, MPI_DOUBLE, allcomptimemultTrans, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&computeTimeMult_, 1, MPI_DOUBLE, allcomptimemult, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (sg::parallel::myGlobalMPIComm->getMyRank() == 0) {
        for (size_t i = 0; i < sg::parallel::myGlobalMPIComm->getNumRanks(); i++) {
          std::cout << "waiting time for process " << i << " is: " << allwaittime[i] << std::endl;
        }

        for (size_t i = 0; i < sg::parallel::myGlobalMPIComm->getNumRanks(); i++) {
          std::cout << "comp time multtrans for process " << i << " is: " << allcomptimemultTrans[i] << std::endl;
        }

        for (size_t i = 0; i < sg::parallel::myGlobalMPIComm->getNumRanks(); i++) {
          std::cout << "comp time mult for process " << i << " is: " << allcomptimemult[i] << std::endl;
        }
      }

      delete[] allwaittime;
      delete[] allcomptimemultTrans;
      delete[] allcomptimemult;
    }

    void DMSystemMatrixVectorizedIdentityMPI::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
#ifdef X86_MIC_SYMMETRIC
      myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);
#endif
      sg::base::DataVector temp(this->numPatchedTrainingInstances_);

      // Operation B
      multVec(alpha, temp);

      // patch result -> set additional entries zero
      if (this->numTrainingInstances_ != temp.getSize()) {
        for (size_t i = 0; i < (temp.getSize() - this->numTrainingInstances_); i++) {
          temp.set(temp.getSize() - (i + 1), 0.0f);
        }
      }

      multTransposeVec(temp, result);

      result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
    }

    void DMSystemMatrixVectorizedIdentityMPI::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
      sg::base::DataVector myClasses(classes);

      // Apply padding
      if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
        myClasses.resizeZero(this->numPatchedTrainingInstances_);
      }

      multTransposeVec(myClasses, b);
    }

    void DMSystemMatrixVectorizedIdentityMPI::rebuildLevelAndIndex() {
      sg::parallel::PartitioningTool::calcDistribution(m_grid.getStorage()->size(), sg::parallel::myGlobalMPIComm->getNumRanks(), _mpi_grid_sizes, _mpi_grid_offsets, 1);
      size_t mpi_rank = sg::parallel::myGlobalMPIComm->getMyRank();

      this->B_->rebuildLevelAndIndex(_mpi_grid_offsets[mpi_rank],
                                     _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank]);
    }

    void DMSystemMatrixVectorizedIdentityMPI::multVec(base::DataVector& alpha, base::DataVector& result) {
      this->myTimer_->start();

      //double funComputationTime = this->B_->multVectorized(alpha, result);
      this->B_->multVectorized(alpha, result);
      //this->computeTimeMult_ += funComputationTime;

      myGlobalMPIComm->Barrier();
      this->computeTimeMult_ += this->myTimer_->stop();
      //double waitTime = this->myTimer_->stop() - computationTime;
      sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_data_offsets, _mpi_data_sizes);

      double completeTime = this->myTimer_->stop();
      this->completeTimeMult_ += completeTime;
      //  double communicationTime = completeTime - computationTime;

      //  double maxComputationTime, minComputationTime;
      //  double maxCompleteTime, minCompleteTime;
      //  double maxCommunicationTime, minCommunicationTime;

      //  MPI_Reduce(&computationTime, &maxComputationTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      //  MPI_Reduce(&computationTime, &minComputationTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      //  MPI_Reduce(&communicationTime, &maxCommunicationTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      //  MPI_Reduce(&communicationTime, &minCommunicationTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      //  MPI_Reduce(&completeTime, &maxCompleteTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      //  MPI_Reduce(&completeTime, &minCompleteTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      //  if (sg::parallel::myGlobalMPIComm->getMyRank() == 0 && false) {
      //    std::cout << "computation     time min - max: " << minComputationTime << " - " << maxComputationTime << " (difference: " << (maxComputationTime - minComputationTime) << ") " << std::endl;
      //    std::cout << "complete        time min - max: " << minCompleteTime << " - " << maxCompleteTime << " (difference: " << (maxCompleteTime - minCompleteTime) << ") " << std::endl;
      //    std::cout << "communication   time min - max: " << minCommunicationTime << " - " << maxCommunicationTime << " (difference: " << (maxCommunicationTime - minCommunicationTime) << ") " << std::endl;
      //    std::cout << std::endl;
      //  }
    }

    void DMSystemMatrixVectorizedIdentityMPI::multTransposeVec(base::DataVector& source, base::DataVector& result) {
      myGlobalMPIComm->Barrier();
      this->myTimer_->start();
      //this->computeTimeMultTrans_ +=
      this->B_->multTransposeVectorized(source, result);

      myGlobalMPIComm->Barrier();
      double comp_time = this->myTimer_->stop();
      this->computeTimeMultTrans_ += comp_time;
      double comp_wait_time = this->myTimer_->stop();
      waitting_time += comp_wait_time - comp_time;

      sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_grid_offsets, _mpi_grid_sizes);

      this->completeTimeMultTrans_ += this->myTimer_->stop();
    }

  }
}
