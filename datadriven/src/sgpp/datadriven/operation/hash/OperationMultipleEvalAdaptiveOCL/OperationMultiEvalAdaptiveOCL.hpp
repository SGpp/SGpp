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
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalAdaptiveOCL/AdaptiveOCLKernelImpl.hpp>

namespace SGPP {
namespace datadriven {

template<typename T>
class OperationMultiEvalAdaptiveOCL: public base::OperationMultipleEval {
 protected:
  size_t m_dims;
  SGPP::base::DataMatrix preparedDataset;
  base::OCLOperationConfiguration parameters;
  T* kernelDataset = nullptr;
  size_t datasetSize = 0;

  T* alphas = nullptr;
  size_t gridSize = 0;
  size_t numSubspaces = 0;

  uint32_t* m_metaInfo;
  T* m_streamingArray;
  T* m_subspaceArray;
  size_t m_streamingCounter;
  size_t m_subspaceCounter;
  size_t m_metaCounter = 0;
  uint32_t m_streamElementCount = 0;
  uint32_t m_subElementCount = 0;

  float m_softAdaptivityLimit = 1;
  uint32_t m_hardAdaptivityLimit = 1;

  /// Timer object to handle time measurements
  SGPP::base::SGppStopwatch myTimer;

  base::GridStorage* storage;

  float_t duration;

  std::shared_ptr<base::OCLManager> manager;
  std::unique_ptr<AdaptiveOCLKernelImpl<T>> kernel;
 public:

  OperationMultiEvalAdaptiveOCL(base::Grid& grid, base::DataMatrix& dataset,
                                std::shared_ptr<base::OCLOperationConfiguration> parameters) :
    OperationMultipleEval(grid, dataset), preparedDataset(dataset),
    parameters(*parameters), myTimer(
      SGPP::base::SGppStopwatch()), duration(-1.0) {
    this->manager = std::make_shared<base::OCLManager>(parameters);

    this->m_dims = dataset.getNcols(); //be aware of transpose!
    this->kernel = std::unique_ptr<AdaptiveOCLKernelImpl<T>>
                   (new AdaptiveOCLKernelImpl<T>(m_dims, this->manager, parameters));

    this->storage = grid.getStorage();
    this->gridSize = grid.getSize();
    this->padDataset(this->preparedDataset);
    this->preparedDataset.transpose();
    this->datasetSize = this->preparedDataset.getNcols();

    //TODO: iteration one -> not a subspace level list, but instead a (level -> points to start of first surplus in surplus array)

    //    std::cout << "dims: " << this->dims << std::endl;
    //    std::cout << "padded instances: " << this->datasetSize << std::endl;

    this->kernelDataset = new T[this->preparedDataset.getNrows() *
                                this->preparedDataset.getNcols()];

    for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
      this->kernelDataset[i] = (T) this->preparedDataset[i];
    }

    m_softAdaptivityLimit = (float)(
                              this->parameters["ADAPTIVE_STREAMING_DENSITY"].getUInt()) / 100.f;
    m_hardAdaptivityLimit = static_cast<uint32_t>
                            (this->parameters["ADAPTIVE_STREAMING_HARD_LIMIT"].getUInt());

    //create the kernel specific data structures
    this->prepare();
  }

  ~OperationMultiEvalAdaptiveOCL() {


    if (m_streamingArray != nullptr) {
      delete m_streamingArray;
    }

    if (m_subspaceArray != nullptr) {
      delete m_subspaceArray;
    }

    if (m_metaInfo != nullptr) {
      delete m_metaInfo;
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


    T* alphaArray = new T[this->gridSize];

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T) alpha[i];
      //printf("i: %i alpha: %f \n", i, alpha[i]);
    }

    for (size_t i = alpha.getSize(); i < this->gridSize; i++) {
      alphaArray[i] = 0.0;
    }

    this->alphas = alphaArray;

    T* resultArray = new T[this->datasetSize];

    for (size_t i = 0; i < this->datasetSize; i++) {
      resultArray[i] = 0.0;
    }

    this->kernel->mult(m_streamingArray, m_subspaceArray, m_metaInfo,
                       this->gridSize, this->kernelDataset, this->datasetSize,
                       this->alphas, resultArray, gridFrom, gridTo, datasetFrom, datasetTo,
                       this->numSubspaces, m_streamElementCount, m_subElementCount);

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

    T* sourceArray = new T[this->datasetSize];

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

    T* resultArray = new T[this->gridSize];

    for (size_t i = 0; i < this->gridSize; i++) {
      resultArray[i] = 0.0;
    }

    this->kernel->multTranspose(m_streamingArray, m_subspaceArray, m_metaInfo,
                                this->gridSize, this->kernelDataset,
                                this->preparedDataset.getNcols(), this->alphas, resultArray, gridFrom, gridTo,
                                datasetFrom, datasetTo,
                                this->numSubspaces, m_streamElementCount, m_subElementCount);

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
    this->buildDatastructure();

    this->kernel->resetKernel();

    //    std::cout << "gridSize: " << this->gridSize << std::endl;
  }

 private:

  size_t padDataset(
    SGPP::base::DataMatrix& dataset) {

    size_t dataBlocking = parameters["LOCAL_SIZE"].getUInt();

    size_t vecWidth = dataBlocking;

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

  uint32_t getSubspaceSize ( uint32_t* level) {

    uint32_t sum = 0;

    //somOf(2^level[d]/2)
    for (size_t d = 0; d < m_dims; d++) {
      if ( level[d] < 1 )
        throw std::runtime_error("Negative level value.");

      //size += static_cast<uint32_t>(powf(2.0f,(float)level[d])/2.0f);
      sum += level[d] - 1;
    }

    return static_cast<uint32_t>(pow(2.0, sum));
  }

  uint32_t calcLinearIndex ( uint32_t* level, uint32_t* index) {
    //calc linear index of gridpoint (relative to subspace)
    uint32_t result = 0;
    uint32_t index_half = 0;
    uint32_t level_calc = 0;

    for (size_t d = 0; d < m_dims - 1; d++) {
      index_half = index[d] >> 1;
      level_calc = static_cast<uint32_t>(pow(2.0, level[d + 1] - 1.0));
      result += index_half;
      result *= level_calc;
    }

    result += index[m_dims - 1] >> 1;

    return result;
  }

  bool isEvalTypeStreaming (uint32_t* level, uint32_t numIndices) {
    uint32_t subSize = getSubspaceSize(level);

    if ( subSize <= m_hardAdaptivityLimit ) {
      return true;
    } else {
      float subspaceDensity = ((float)numIndices) / ((float)subSize);

      if (subspaceDensity <= m_softAdaptivityLimit) {
        return true;
      } else {
        return false;
      }
    }
  }

  void addMetaEntry(uint32_t* level, uint32_t startIndex, uint32_t numIndices,
                    bool isStreaming) {
    size_t metaInfoStepSize = m_dims + 3;
    m_metaInfo[metaInfoStepSize * m_metaCounter] = isStreaming ? 1 : 0;
    m_metaInfo[metaInfoStepSize * m_metaCounter + 1] = startIndex;
    m_metaInfo[metaInfoStepSize * m_metaCounter + 2] = numIndices;

    for ( int d = 0; d < m_dims; d++) {
      m_metaInfo[metaInfoStepSize * m_metaCounter + 3 + d] = level[d];
    }

    m_metaCounter++;
  }

  void addToStreamingArray(uint32_t* level, std::vector<uint32_t> indexList,
                           uint32_t numIndices) {
    //read and store indices
    for (int i = 0; i < numIndices; i++) {

      size_t indexStepSize = (m_dims + 1); //indexVector + gridIndex

      //store index
      for (int d = 0; d < m_dims; d++) {
        T tmpVal = (T) (indexList[(m_dims) + (i * (m_dims + 1)) + d + 1]);
        m_streamingArray[(m_streamingCounter + i)*indexStepSize + d] = tmpVal;
      }

      //store gridIndex
      T gridIndex = (T) (indexList[(this->m_dims) + (i * (this->m_dims + 1))]);
      m_streamingArray[(m_streamingCounter + i)*indexStepSize + m_dims] = gridIndex;
    }

    addMetaEntry(level, static_cast<uint32_t>(m_streamingCounter), numIndices,
                 true);
    m_streamingCounter += numIndices;
  }

  void addToSubspaceArray(uint32_t* level, std::vector<uint32_t> indexList,
                          uint32_t numIndices) {
    uint32_t subSize = getSubspaceSize(level);
    uint32_t indexStepSize = static_cast<uint32_t>(m_dims + 1);

    for (int i = 0; i < subSize; i++) {
      for (int d = 0; d < m_dims; d++) {
        m_subspaceArray[(m_subspaceCounter + i)*indexStepSize + d] = static_cast<T>(0);
      }

      m_subspaceArray[(m_subspaceCounter + i)*indexStepSize + m_dims] =
        static_cast<T>(NAN);
    }

    //read and store indices
    for (int i = 0; i < numIndices; i++) {

      uint32_t* index = new uint32_t[m_dims];

      //store index
      int levelIndexSum = 0;

      for (int d = 0; d < m_dims; d++) {
        T tmpVal = (T) (indexList[(m_dims) + (i * (m_dims + 1)) + d + 1]);
        index[d] = static_cast<uint32_t>(tmpVal);
        levelIndexSum += index[d] + level[d];
      }

      size_t linIndex = calcLinearIndex(level, index);

      for (int d = 0; d < m_dims; d++) {
        m_subspaceArray[(m_subspaceCounter + linIndex)*indexStepSize + d] =
          static_cast<T>(index[d]);
      }

      //alternative
      //memcpy(&m_subspaceArray[(m_subspaceCounter + linIndex)*indexStepSize], index, sizeof(uint32_t)*m_dims);

      //store gridIndex
      T gridIndex = (T) (indexList[(this->m_dims) + (i * (this->m_dims + 1))]);

      //ugly hack to fix padding points messing up the calculation
      //TODO: REMOVE!!
      if ( levelIndexSum == m_dims * 2) {
        gridIndex = 0;
      }

      m_subspaceArray[(m_subspaceCounter + linIndex)*indexStepSize + m_dims] =
        gridIndex;
    }

    addMetaEntry(level, static_cast<uint32_t>(m_subspaceCounter), subSize, false);
    m_subspaceCounter += subSize;
  }

  void buildDatastructure() {
    size_t dataBlocking = parameters["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
    size_t transGridBlocking =
      parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"].getUInt();

    //TODO: is this a bug, Raphael? (David)
    size_t blockingSize = std::max(dataBlocking, transGridBlocking);

    uint32_t localWorkSize = static_cast<uint32_t>
                             (parameters["LOCAL_SIZE"].getUInt()) * static_cast<uint32_t>(blockingSize);

    size_t remainder = this->storage->size() % localWorkSize;
    size_t padding = 0;

    if (remainder != 0) {
      padding = localWorkSize - remainder;
    }

    this->gridSize = this->storage->size() + padding;

    m_streamingCounter = 0;
    m_subspaceCounter = 0;
    m_metaCounter = 0;
    m_streamElementCount = 0;
    m_subElementCount = 0;

    uint32_t maxLevel = (uint32_t) (this->storage->getMaxLevel());

    std::map<uint32_t, std::vector<uint32_t> > flatLevelList;
    flatLevelList.clear();

    this->numSubspaces = 0;

    //vector used to preserve subspace order as std::map silently re-arranges items
    std::vector<uint32_t> flatLevelOrder;

    uint32_t curLevel = 1;
    uint32_t curIndex = 1;

    //iterate through the grid to read out levels and indices
    for (size_t gridIndex = 0; gridIndex < this->gridSize; gridIndex++) {
      uint32_t* level = new uint32_t[this->m_dims];
      uint32_t* index = new uint32_t[this->m_dims];

      if (gridIndex < this->storage->size()) {
        SGPP::base::GridIndex* point = this->storage->get(gridIndex);

        for (size_t d = 0; d < this->m_dims; d++) {
          point->get(d, curLevel, curIndex);
          level[d] = curLevel;
          index[d] = curIndex;
        }
      }
      //add points for padding that don't affect evaluation
      else {
        for (size_t d = 0; d < this->m_dims; d++) {
          level[d] = 1;
          index[d] = 1;
        }
      }


      //calculate flatLevel to uniquely identify each subspace
      uint32_t flatLevel = flattenLevel(this->m_dims, level, maxLevel);
      std::map<uint32_t, std::vector<uint32_t> >::iterator flatLevelEntry =
        flatLevelList.find(flatLevel);

      //flatLevel doesn't exist, add a new subspace
      if (flatLevelEntry == flatLevelList.end()) {
        std::vector<uint32_t> indexList;
        flatLevelOrder.push_back(flatLevel);

        for (size_t d = 0; d < m_dims; d++) {
          indexList.push_back(level[d]);
        }

        /*printf("new subspace added! flatLevel: %i \n", flatLevel);
         printf("level: %i \t %i \t %i \t %i \n", level[0], level[1], level[2], level[3]);
         printf("index: %i \t %i \t %i \t %i \n", index[0], index[1], index[2], index[3]);
         printf("gridindex: %i \n", gridIndex);*/

        //store gridIndex to map to alphaArray
        indexList.push_back((uint32_t) gridIndex);

        //store index
        for (size_t d = 0; d < this->m_dims; d++) {
          indexList.push_back(index[d]);
        }

        flatLevelList.insert(std::make_pair(flatLevel, indexList));
        this->numSubspaces += 1;
      } else {
        //store gridIndex to map to alphaArray
        flatLevelEntry->second.push_back((uint32_t) gridIndex);

        //store index
        for (size_t d = 0; d < this->m_dims; d++) {
          flatLevelEntry->second.push_back(index[d]);
        }

        /*printf("level: %i \t %i \t %i \t %i \n", flatLevelEntry->second[0], flatLevelEntry->second[1], flatLevelEntry->second[2], flatLevelEntry->second[3]);
         printf("index: %i \t %i \t %i \t %i \n", index[0], index[1], index[2], index[3]);
         printf("gridindex: %i \n", gridIndex);*/
      }
    }



    //    flatLevelList.clear();

    printf("GridPoints: %lu Subspaces: %lu \n", this->gridSize, numSubspaces);
    printf("Subspace Utilization: %f \n",
           ((double)(this->gridSize)) / ((double)numSubspaces));


    //TODO: dont repeat the loop
    //iterate our previously constructed contains and fill our arrays
    typedef std::map<uint32_t, std::vector<uint32_t> >::iterator it_type;

    for (auto& flatLevel : flatLevelOrder) {
      //as mentioned above this preserves the initial order of appearance
      it_type kv = flatLevelList.find(flatLevel);

      uint32_t numIndices = static_cast<uint32_t>((static_cast<size_t>
                            (kv->second.size()) - this->m_dims) / (this->m_dims + 1));

      //read the level
      uint32_t* level = new uint32_t[this->m_dims];

      for (int d = 0; d < this->m_dims; d++) {
        level[d] = kv->second[d];
      }

      if ( isEvalTypeStreaming(level, numIndices) ) {

        m_streamElementCount += numIndices;
      } else {
        m_subElementCount += getSubspaceSize(level);
      }

    }

    uint32_t metaInfoSize = static_cast<uint32_t>(m_dims) +
                            3; //level vector + start, size and type
    uint32_t subElementSize = static_cast<uint32_t>(m_dims) + 1;
    uint32_t streamElementSize = static_cast<uint32_t>(m_dims) + 1;

    this->m_metaInfo = new uint32_t[this->numSubspaces * metaInfoSize];
    this->m_subspaceArray = new T[m_subElementCount * subElementSize];
    this->m_streamingArray = new T[m_streamElementCount * streamElementSize];

    for (auto& flatLevel : flatLevelOrder) {
      //as mentioned above this preserves the initial order of appearance
      it_type kv = flatLevelList.find(flatLevel);

      uint32_t numIndices = static_cast<uint32_t>((kv->second.size() -
                            this->m_dims) / (this->m_dims + 1));

      //read the level
      uint32_t* level = new uint32_t[this->m_dims];

      for (int d = 0; d < this->m_dims; d++) {
        level[d] = kv->second[d];
      }

      if ( isEvalTypeStreaming(level, numIndices) ) {

        addToStreamingArray(level, kv->second, numIndices);
      } else {
        addToSubspaceArray(level, kv->second, numIndices);
      }

    }

    printf("subElementCount: %i, streamElementCount: %i \n", m_subElementCount,
           m_streamElementCount);


  }
};
}
}
