// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPI.hpp>
#include <sgpp/parallel/operation/SPParallelOpFactory.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

DMSystemMatrixSPVectorizedIdentityMPI::DMSystemMatrixSPVectorizedIdentityMPI(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataMatrixSP& trainData, float lambda,
    VectorizationType vecMode)
    : DMSystemMatrixBaseSP(trainData, lambda),
      vecMode_(vecMode),
      numTrainingInstances_(0),
      numPatchedTrainingInstances_(0),
      m_grid(SparseGrid) {
  // handle unsupported vector extensions
  if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC &&
      this->vecMode_ != OpenCL &&
      this->vecMode_ != Hybrid_X86SIMD_OpenCL) {
    throw sgpp::base::operation_exception(
        "DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
  }

  this->dataset_ = new sgpp::base::DataMatrixSP(trainData);
  this->numTrainingInstances_ = this->dataset_->getNrows();
  this->numPatchedTrainingInstances_ =
      sgpp::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);

  this->dataset_->transpose();

  size_t mpi_size = sgpp::parallel::myGlobalMPIComm->getNumRanks();
  size_t mpi_rank = sgpp::parallel::myGlobalMPIComm->getMyRank();

  // arrays for distribution settings
  _mpi_grid_sizes = new int[mpi_size];
  _mpi_grid_offsets = new int[mpi_size];
  _mpi_data_sizes = new int[mpi_size];
  _mpi_data_offsets = new int[mpi_size];

  // calculate distribution
  sgpp::parallel::PartitioningTool::calcDistribution(m_grid.getStorage().getSize(), mpi_size,
                                                     _mpi_grid_sizes, _mpi_grid_offsets, 1);
  sgpp::parallel::PartitioningTool::calcDistribution(
      this->numPatchedTrainingInstances_, mpi_size, _mpi_data_sizes, _mpi_data_offsets,
      sgpp::parallel::DMVectorizationPaddingAssistant::getVecWidthSP(this->vecMode_));

  this->B_ = sgpp::op_factory::createOperationMultipleEvalVectorizedSP(
      m_grid, this->vecMode_, this->dataset_, _mpi_grid_offsets[mpi_rank],
      _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank], _mpi_data_offsets[mpi_rank],
      _mpi_data_offsets[mpi_rank] + _mpi_data_sizes[mpi_rank]);
  waitting_time = 0.0;
}

DMSystemMatrixSPVectorizedIdentityMPI::~DMSystemMatrixSPVectorizedIdentityMPI() {
  delete this->B_;
  delete this->dataset_;

  delete[] this->_mpi_grid_sizes;
  delete[] this->_mpi_grid_offsets;
  delete[] this->_mpi_data_sizes;
  delete[] this->_mpi_data_offsets;
  double* allwaittime = new double[sgpp::parallel::myGlobalMPIComm->getNumRanks()];
  double* allcomptimemultTrans = new double[sgpp::parallel::myGlobalMPIComm->getNumRanks()];
  double* allcomptimemult = new double[sgpp::parallel::myGlobalMPIComm->getNumRanks()];
  MPI_Gather(&waitting_time, 1, MPI_DOUBLE, allwaittime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&computeTimeMultTrans_, 1, MPI_DOUBLE, allcomptimemultTrans, 1, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&computeTimeMult_, 1, MPI_DOUBLE, allcomptimemult, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
    for (size_t i = 0; i < sgpp::parallel::myGlobalMPIComm->getNumRanks(); i++) {
      std::cout << "waiting time for process " << i << " is: " << allwaittime[i] << std::endl;
    }

    for (size_t i = 0; i < sgpp::parallel::myGlobalMPIComm->getNumRanks(); i++) {
      std::cout << "comp time multtrans for process " << i << " is: " << allcomptimemultTrans[i]
                << std::endl;
    }

    for (size_t i = 0; i < sgpp::parallel::myGlobalMPIComm->getNumRanks(); i++) {
      std::cout << "comp time mult for process " << i << " is: " << allcomptimemult[i] << std::endl;
    }
  }

  delete[] allwaittime;
  delete[] allcomptimemultTrans;
  delete[] allcomptimemult;
}

void DMSystemMatrixSPVectorizedIdentityMPI::mult(sgpp::base::DataVectorSP& alpha,
                                                 sgpp::base::DataVectorSP& result) {
#ifdef X86_MIC_SYMMETRIC
  myGlobalMPIComm->broadcastSPGridCoefficientsFromRank0(alpha);
#endif
  sgpp::base::DataVectorSP temp(this->numPatchedTrainingInstances_);

  // Operation B
  multVec(alpha, temp);

  // patch result -> set additional entries zero
  if (this->numTrainingInstances_ != temp.getSize()) {
    for (size_t i = 0; i < (temp.getSize() - this->numTrainingInstances_); i++) {
      temp.set(temp.getSize() - (i + 1), 0.0f);
    }
  }

  multTransposeVec(temp, result);

  result.axpy(static_cast<float>(this->numTrainingInstances_) * this->lambda_, alpha);
}

void DMSystemMatrixSPVectorizedIdentityMPI::generateb(sgpp::base::DataVectorSP& classes,
                                                      sgpp::base::DataVectorSP& b) {
  sgpp::base::DataVectorSP myClasses(classes);

  // Apply padding
  if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
    myClasses.resizeZero(this->numPatchedTrainingInstances_);
  }

  multTransposeVec(myClasses, b);
}

void DMSystemMatrixSPVectorizedIdentityMPI::rebuildLevelAndIndex() {
  sgpp::parallel::PartitioningTool::calcDistribution(m_grid.getStorage().getSize(),
                                                     sgpp::parallel::myGlobalMPIComm->getNumRanks(),
                                                     _mpi_grid_sizes, _mpi_grid_offsets, 1);
  size_t mpi_rank = sgpp::parallel::myGlobalMPIComm->getMyRank();

  this->B_->rebuildLevelAndIndex(_mpi_grid_offsets[mpi_rank],
                                 _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank]);
}

void DMSystemMatrixSPVectorizedIdentityMPI::multVec(base::DataVectorSP& alpha,
                                                    base::DataVectorSP& result) {
  this->myTimer_->start();

  this->computeTimeMult_ += this->B_->multVectorized(alpha, result);

  sgpp::parallel::myGlobalMPIComm->dataVectorSPAllToAll(result, _mpi_data_offsets, _mpi_data_sizes);

  this->completeTimeMult_ += this->myTimer_->stop();
}

void DMSystemMatrixSPVectorizedIdentityMPI::multTransposeVec(base::DataVectorSP& source,
                                                             base::DataVectorSP& result) {
  this->myTimer_->start();
  this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(source, result);

  sgpp::parallel::myGlobalMPIComm->dataVectorSPAllToAll(result, _mpi_grid_offsets, _mpi_grid_sizes);

  this->completeTimeMultTrans_ += this->myTimer_->stop();
}

}  // namespace parallel
}  // namespace sgpp
