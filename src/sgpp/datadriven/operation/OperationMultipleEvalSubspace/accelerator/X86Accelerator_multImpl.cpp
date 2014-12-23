#include "../../OperationMultipleEvalSubspace/accelerator/X86Accelerator.hpp"

void X86Accelerator::multImpl(sg::base::DataVector &alpha, sg::base::DataVector &result, const size_t start_index_data,
        const size_t end_index_data) {

    size_t tid = omp_get_thread_num();
    if (tid == 0) {
        this->setCoefficients(result);

        //build bitwise copyable structures

        //create a subspace representation for the MIC
        uint32_t hInverseOffset = 0;
        uint32_t surplusOffset = 0;

        for (X86AcceleratorSubspaceNode &subspace : this->allSubspaceNodes) {
            X86AcceleratorSubspaceNodeMIC micSubspace;
            micSubspace.arriveDiff = subspace.arriveDiff;
            micSubspace.existingGridPointsOnLevel = subspace.existingGridPointsOnLevel;
            micSubspace.flatLevel = subspace.flatLevel;
            micSubspace.hInverseOffset = hInverseOffset;
            for (uint32_t i = 0; i < this->dim; i++) {
                this->allHInverseMic.push_back(subspace.hInverse[i]);
            }
            hInverseOffset += this->dim;
            micSubspace.jumpTargetIndex = subspace.jumpTargetIndex;
            micSubspace.surplusOffset = surplusOffset;
            micSubspace.gridPointsOnLevel = subspace.gridPointsOnLevel;
            // if (subspace.type == X86AcceleratorSubspaceNode::LIST) {
            //     for (uint32_t gridPointIndex = 0; gridPointIndex < subspace.gridPointsOnLevel; gridPointIndex++) {
            //         this->allSurplussesMic.push_back(std::numeric_limits<double>::quiet_NaN());
            //     }
            //     for (pair<uint32_t, double> pair : subspace.indexFlatSurplusPairs) {
            //         this->allSurplussesMic[surplusOffset + pair.first] = pair.second;
            //     }
            // } else {
            // ARRAY representation
            for (uint32_t gridPointIndex = 0; gridPointIndex < subspace.gridPointsOnLevel; gridPointIndex++) {
                //cout << "read: " << subspace.subspaceArray[gridPointIndex] << endl;
                this->allSurplussesMic.push_back(subspace.subspaceArray[gridPointIndex]);

            }
            //}
            surplusOffset += subspace.gridPointsOnLevel;
            //cout << "surplusOffset: " << surplusOffset << endl;
            //cout << "gridPointsOnLevel: " << subspace.gridPointsOnLevel << endl;
            this->allSubspaceNodesMic.push_back(micSubspace);
        }
        //cout << "subspace: " << this->allSubspaceNodes.size() << endl;

    }

#pragma omp barrier

    size_t dim = this->dataset->getNcols();

    const double * const datasetPtr = this->dataset->getPointer();

    size_t dataElements = this->dataset->getNrows() * this->dataset->getNcols();
    size_t totalSurplusCount = this->allSurplussesMic.size();
    double *allSurplussesMicArray = this->allSurplussesMic.data();
    size_t totalHInverseCount = this->allHInverseMic.size();
    uint32_t *allHInverseMicArray = this->allHInverseMic.data();
    X86AcceleratorSubspaceNodeMIC *allSubspaceNodesMicArray = this->allSubspaceNodesMic.data();

//#pragma offload target(mic)	in(datasetPtr:length(dataElements))	inout(allSurplussesMicArray:length(totalSurplusCount)) in(allHInverseMicArray:length(totalHInverseCount) in(allSubspaceNodesMicArray:length(this->subspaceCount)) in(start_index_data, end_index_data)
    {

        size_t totalThreadNumber = X86ACCELERATOR_PARALLEL_DATA_POINTS + X86ACCELERATOR_VEC_PADDING;

        double *evalIndexValuesAll = new double[(dim + 1) * totalThreadNumber];
        for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
            evalIndexValuesAll[i] = 1.0;
        }

        //for faster index flattening
        uint32_t *intermediatesAll = new uint32_t[(dim + 1) * totalThreadNumber];
        for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
            intermediatesAll[i] = 0.0;
        }

        size_t validIndices[X86ACCELERATOR_PARALLEL_DATA_POINTS + X86ACCELERATOR_VEC_PADDING];
        size_t validIndicesCount;

        size_t levelIndices[X86ACCELERATOR_PARALLEL_DATA_POINTS + X86ACCELERATOR_VEC_PADDING];

        //#pragma omp parallel for
        for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; dataIndexBase +=
        X86ACCELERATOR_PARALLEL_DATA_POINTS) {

            for (size_t i = 0; i < totalThreadNumber; i++) {
                levelIndices[i] = 0;
                //nextIterationToRecalcReferences[i] = 0;
            }

            for (size_t subspaceIndex = 0; subspaceIndex < subspaceCount; subspaceIndex++) {
                //X86AcceleratorSubspaceNodeMIC &subspace = this->allSubspaceNodesMic[subspaceIndex];
                X86AcceleratorSubspaceNodeMIC &subspace = allSubspaceNodesMicArray[subspaceIndex];

                validIndicesCount = 0;
                for (size_t parallelIndex = 0; parallelIndex < X86ACCELERATOR_PARALLEL_DATA_POINTS; parallelIndex++) {
                    size_t parallelLevelIndex = levelIndices[parallelIndex];

                    if (parallelLevelIndex == subspaceIndex) {
                        validIndices[validIndicesCount] = parallelIndex;
                        validIndicesCount += 1;
                    }
                }

                //padding for up to vector size, no padding required if all data tuples participate as
                //the number of data points is a multiple of the vector width
                size_t paddingSize = min((int) (validIndicesCount + X86ACCELERATOR_VEC_PADDING),
                X86ACCELERATOR_PARALLEL_DATA_POINTS + X86ACCELERATOR_VEC_PADDING);
                for (size_t i = validIndicesCount; i < paddingSize; i++) {
                    size_t threadId = X86ACCELERATOR_PARALLEL_DATA_POINTS + (i - validIndicesCount);
                    validIndices[i] = threadId;
                    levelIndices[threadId] = 0;
                    //nextIterationToRecalcReferences[threadId] = 0;
                    double *evalIndexValues = evalIndexValuesAll + (dim + 1) * threadId;

                    //for faster index flattening, last element is for padding
                    uint32_t *intermediates = intermediatesAll + (dim + 1) * threadId;
                    for (size_t j = 0; j < dim; j++) {
                        evalIndexValues[j] = 1.0;
                        intermediates[j] = 0;
                    }
                }

                listMultInner(dim, datasetPtr, alpha, dataIndexBase, end_index_data, subspace,
                        &(allSurplussesMicArray[subspace.surplusOffset]), validIndicesCount, validIndices, levelIndices,
                        evalIndexValuesAll, intermediatesAll, &(allHInverseMicArray[subspace.hInverseOffset]));

            } // end iterate subspaces
        } // end iterate chunks
        delete evalIndexValuesAll;
        delete intermediatesAll;
    }

#pragma omp barrier

    if (tid == 0) {
        //copy the accelerator data structure back into the host structure

        for (size_t i = 0; i < this->allSubspaceNodes.size(); i++) {
            //TODO: currently only considers ARRAY representation
            X86AcceleratorSubspaceNode &subspace = this->allSubspaceNodes[i];
            X86AcceleratorSubspaceNodeMIC &micSubspace = this->allSubspaceNodesMic[i];
            // ARRAY representation
            // cout << "gridPointsOnLevel: " << subspace.gridPointsOnLevel << endl;
            for (uint32_t gridPointIndex = 0; gridPointIndex < subspace.gridPointsOnLevel; gridPointIndex++) {
                // cout << "gridPointIndex: " << gridPointIndex << endl;
                // cout << "value: " << this->allSurplussesMic[micSubspace.surplusOffset + gridPointIndex] << endl;
                // cout << "value reference: " << subspace.subspaceArray[gridPointIndex] << endl;
                //this->allSurplussesMic.push_back(subspace.subspaceArray[gridPointIndex]);
                subspace.subspaceArray[gridPointIndex] = this->allSurplussesMic[micSubspace.surplusOffset
                        + gridPointIndex];
            }
        }

        this->unflatten(result);
        this->allSurplussesMic.clear();
        this->allHInverseMic.clear();
        this->allSubspaceNodesMic.clear();
    }
}
