// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include "../../../opencl/OCLConfigurationParameters.hpp"
#include "../../../opencl/OCLManager.hpp"
#include "AdaptiveOCLKernelImpl.hpp"

namespace SGPP {
namespace datadriven {

template<typename T>
class OperationMultiEvalAdaptiveOCL: public base::OperationMultipleEval {
protected:
    size_t dims;SGPP::base::DataMatrix preparedDataset;
    base::OCLConfigurationParameters parameters;
    T *kernelDataset = nullptr;
    size_t datasetSize = 0;
    /// Member to store the sparse grid's subspaces
    //SubspaceInfo* subSpaces = nullptr;
    std::map<uint32_t, std::vector<uint32_t>> subspaceMap;

    T* levels = nullptr;
    T* indices = nullptr;
    T* alphas = nullptr;
    /// Member to store the sparse grid's indices for better vectorization
    T* index = nullptr;
    size_t gridSize = 0;
    size_t numSubspaces = 0;
    /// Timer object to handle time measurements
    SGPP::base::SGppStopwatch myTimer;

    base::GridStorage* storage;

    float_t duration;

    std::shared_ptr<base::OCLManager> manager;
    std::unique_ptr<AdaptiveOCLKernelImpl<T>> kernel;
public:

    OperationMultiEvalAdaptiveOCL(base::Grid& grid, base::DataMatrix& dataset,
            std::shared_ptr<base::OCLConfigurationParameters> parameters) :
            OperationMultipleEval(grid, dataset), preparedDataset(dataset), parameters(*parameters), myTimer(
            SGPP::base::SGppStopwatch()), duration(-1.0) {
        this->manager = std::make_shared<base::OCLManager>(parameters);

        this->dims = dataset.getNcols(); //be aware of transpose!
        this->kernel = std::unique_ptr<AdaptiveOCLKernelImpl<T>>(new AdaptiveOCLKernelImpl<T>(dims, this->manager, parameters));

        this->storage = grid.getStorage();
        this->gridSize = grid.getSize();
        this->padDataset(this->preparedDataset);
        this->preparedDataset.transpose();
        this->datasetSize = this->preparedDataset.getNcols();

        //TODO: iteration one -> not a subspace level list, but instead a (level -> points to start of first surplus in surplus array)

//    std::cout << "dims: " << this->dims << std::endl;
//    std::cout << "padded instances: " << this->datasetSize << std::endl;

        this->kernelDataset = new T[this->preparedDataset.getNrows() * this->preparedDataset.getNcols()];
        for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
            this->kernelDataset[i] = (T) this->preparedDataset[i];
        }

        //create the kernel specific data structures
        this->prepare();
    }

    ~OperationMultiEvalAdaptiveOCL() {

        if (this->index != nullptr) {
            delete this->index;
        }

        if (this->kernelDataset != nullptr) {
            delete this->kernelDataset;
        }
    }

    void mult(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result) override {
        this->myTimer.start();

        size_t gridFrom = 0;
        size_t gridTo = this->gridSize;
        size_t datasetFrom = 0;
        size_t datasetTo = this->datasetSize;

        T *alphaArray = new T[this->gridSize];
        for (size_t i = 0; i < alpha.getSize(); i++) {
            alphaArray[i] = (T) alpha[i];
        }
        for (size_t i = alpha.getSize(); i < this->gridSize; i++) {
            alphaArray[i] = 0.0;
        }
        this->alphas = alphaArray;
        //recalculateLevelAndIndex();

        T *resultArray = new T[this->datasetSize];
        for (size_t i = 0; i < this->datasetSize; i++) {
            resultArray[i] = 0.0;
        }

        this->kernel->mult(this->levels, this->indices, this->gridSize, this->kernelDataset, this->datasetSize,
                this->alphas, resultArray, gridFrom, gridTo, datasetFrom, datasetTo, this->numSubspaces);

        for (size_t i = 0; i < result.getSize(); i++) {
            result[i] = resultArray[i];
        }

        delete alphaArray;
        delete resultArray;
        this->duration = this->myTimer.stop();
    }

    void multTranspose(
    SGPP::base::DataVector& source,
    SGPP::base::DataVector& result) override {
        this->myTimer.start();

        size_t gridFrom = 0;
        size_t gridTo = this->gridSize;
        size_t datasetFrom = 0;
        size_t datasetTo = this->datasetSize;

        T *sourceArray = new T[this->datasetSize];
        for (size_t i = 0; i < source.getSize(); i++) {
            sourceArray[i] = (T) source[i];
            //printf("source: %f \n", sourceArray[i]);
        }
        for (size_t i = source.getSize(); i < this->datasetSize; i++) {
            sourceArray[i] = 0.0;
            //printf("source: %f \n", sourceArray[i]);
        }
        this->alphas = sourceArray;
        //recalculateLevelAndIndex();

        T *resultArray = new T[this->gridSize];
        for (size_t i = 0; i < this->gridSize; i++) {
            resultArray[i] = 0.0;
        }

        this->kernel->multTranspose(this->levels, this->indices, this->gridSize, this->kernelDataset,
                this->preparedDataset.getNcols(), this->alphas, resultArray, gridFrom, gridTo, datasetFrom, datasetTo,
                this->numSubspaces);

        for (size_t i = 0; i < result.getSize(); i++) {
            result[i] = resultArray[i];
        }

        delete sourceArray;
        delete resultArray;
        this->duration = this->myTimer.stop();
    }

    float_t getDuration() {
        return this->duration;
    }

    void prepare() override {
        this->recalculateLevelAndIndex();

        this->kernel->resetKernel();

//    std::cout << "gridSize: " << this->gridSize << std::endl;
    }

private:

    size_t padDataset(
    SGPP::base::DataMatrix& dataset) {

        size_t dataBlocking = parameters.getAsUnsigned("KERNEL_DATA_BLOCKING_SIZE");
        size_t transGridBlocking = parameters.getAsUnsigned("KERNEL_TRANS_GRID_BLOCKING_SIZE");

        size_t blockingSize = std::max(dataBlocking, transGridBlocking);

        size_t vecWidth = parameters.getAsUnsigned("LOCAL_SIZE") * blockingSize;

        // Assure that data has a even number of instances -> padding might be needed
        size_t remainder = dataset.getNrows() % vecWidth;
        size_t loopCount = vecWidth - remainder;

        if (loopCount != vecWidth) {
            SGPP::base::DataVector lastRow(dataset.getNcols());
            size_t oldSize = dataset.getNrows();
            dataset.getRow(dataset.getNrows() - 1, lastRow);
            dataset.resize(dataset.getNrows() + loopCount);

            for (size_t i = 0; i < loopCount; i++) {
                dataset.setRow(oldSize + i, lastRow);
            }
        }

        return dataset.getNrows();
    }

    uint32_t flattenLevel(size_t dims, uint32_t* level, size_t maxLevel) {
        uint32_t levelFlat = 0;

        levelFlat += level[dims - 1];

        // loop terminates at -1
        for (int i = static_cast<int>(dims - 2); i >= 0; i--) {
            levelFlat *= static_cast<uint32_t>(maxLevel);
            levelFlat += level[i];
        }

        return levelFlat;
    }

    uint32_t flattenIndex(size_t dim, uint32_t* maxIndices, uint32_t* index) {
        uint32_t indexFlat = index[0];
        indexFlat >>= 1;

        for (size_t i = 1; i < dim; i++) {
            uint32_t actualDirectionGridPoints = maxIndices[i];
            actualDirectionGridPoints >>= 1;
            indexFlat *= actualDirectionGridPoints;
            uint32_t actualIndex = index[i];
            actualIndex >>= 1; //divide index by 2, skip even indices
            indexFlat += actualIndex;
        }

        return indexFlat;
    }

    void recalculateLevelAndIndex() {

        if (this->index != nullptr)
            delete this->index;

        size_t dataBlocking = parameters.getAsUnsigned("KERNEL_DATA_BLOCKING_SIZE");
        size_t transGridBlocking = parameters.getAsUnsigned("KERNEL_TRANS_GRID_BLOCKING_SIZE");

        //TODO: is this a bug, Raphael? (David)
        size_t blockingSize = std::max(dataBlocking, transGridBlocking);

        uint32_t localWorkSize = static_cast<uint32_t>(parameters.getAsUnsigned("LOCAL_SIZE")) * static_cast<uint32_t>(blockingSize);

        size_t remainder = this->storage->size() % localWorkSize;
        size_t padding = 0;

        if (remainder != 0) {
            padding = localWorkSize - remainder;
        }

        this->gridSize = this->storage->size() + padding;

        uint32_t maxLevel = (uint32_t) (this->storage->getMaxLevel());
        uint32_t curLevel = 1;
        uint32_t curIndex = 1;

        std::map<uint32_t, std::vector<uint32_t> > flatLevelList;
        flatLevelList.clear();

        this->numSubspaces = 0;

        std::vector<uint32_t> flatLevelOrder;

        for (size_t gridIndex = 0; gridIndex < this->gridSize; gridIndex++) {
            SGPP::base::GridIndex* point = this->storage->get(gridIndex);

            uint32_t* level = new uint32_t[this->dims];
            uint32_t* index = new uint32_t[this->dims];

            if (gridIndex < this->storage->size()) {
                for (size_t d = 0; d < this->dims; d++) {
                    point->get(d, curLevel, curIndex);
                    level[d] = curLevel;
                    index[d] = curIndex;
                }
            }
            //add points for padding that dont affect evaluation
            else {
                for (size_t d = 0; d < this->dims; d++) {
                    level[d] = 1.0f;
                    index[d] = 1.0f;
                }
            }

            uint32_t flatLevel = flattenLevel(this->dims, level, maxLevel);
            std::map<uint32_t, std::vector<uint32_t> >::iterator flatLevelEntry = flatLevelList.find(flatLevel);

            if (flatLevelEntry == flatLevelList.end()) {
                std::vector<uint32_t> indexList;
                flatLevelOrder.push_back(flatLevel);

                for (size_t d = 0; d < this->dims; d++) {
                    indexList.push_back(level[d]);
                }
                /*printf("new subspace added! flatLevel: %i \n", flatLevel);
                 printf("level: %i \t %i \t %i \t %i \n", level[0], level[1], level[2], level[3]);
                 printf("index: %i \t %i \t %i \t %i \n", index[0], index[1], index[2], index[3]);
                 printf("gridindex: %i \n", gridIndex);*/

                indexList.push_back((uint32_t) gridIndex);

                for (size_t d = 0; d < this->dims; d++) {
                    indexList.push_back(index[d]);
                }

                flatLevelList.insert(std::make_pair(flatLevel, indexList));
                this->numSubspaces += 1;
            } else {
                flatLevelEntry->second.push_back((uint32_t) gridIndex);
                for (size_t d = 0; d < this->dims; d++) {
                    flatLevelEntry->second.push_back(index[d]);
                }

                /*printf("level: %i \t %i \t %i \t %i \n", flatLevelEntry->second[0], flatLevelEntry->second[1], flatLevelEntry->second[2], flatLevelEntry->second[3]);
                 printf("index: %i \t %i \t %i \t %i \n", index[0], index[1], index[2], index[3]);
                 printf("gridindex: %i \n", gridIndex);*/
            }
        }

//    flatLevelList.clear();

        this->levels = new T[(this->dims + 1) * numSubspaces];
        this->indices = new T[(this->dims + 1) * this->gridSize];

        size_t indexCounter = 0;
        size_t levelCounter = 0;

        typedef std::map<uint32_t, std::vector<uint32_t> >::iterator it_type;
        for (auto &flatLevel : flatLevelOrder) {

            it_type kv = flatLevelList.find(flatLevel);

            size_t numIndices = (kv->second.size() - this->dims) / (this->dims + 1);

            //store the number of corresponding indices in front of the level-vector
            this->levels[(levelCounter * (this->dims + 1))] = (T) (indexCounter + numIndices);

            //read and store level
            for (int d = 0; d < this->dims; d++) {

                this->levels[(levelCounter * (this->dims + 1)) + d + 1] = (T) (kv->second[d]);
            }

            //printf("flat_level: %i \n", kv->first);
            //printf("level_stored: %i \t %i \t %i \t %i \n", kv->second[0], kv->second[1], kv->second[2], kv->second[3]);

            //read and store indices
            for (int i = 0; i < numIndices; i++) {

                //store alpha index
                T alpha_idx = (T) (kv->second[(this->dims) + (i * (this->dims + 1))]);
                this->indices[(indexCounter + i) * (this->dims + 1)] = alpha_idx;

                for (int d = 0; d < this->dims; d++) {
                    T tmpVal = (T) (kv->second[(this->dims) + (i * (this->dims + 1)) + d + 1]);
                    this->indices[(indexCounter + i) * (this->dims + 1) + d + 1] = tmpVal;
                }
            }

            indexCounter += numIndices;
            levelCounter += 1;
        }

    }

};

}
}
