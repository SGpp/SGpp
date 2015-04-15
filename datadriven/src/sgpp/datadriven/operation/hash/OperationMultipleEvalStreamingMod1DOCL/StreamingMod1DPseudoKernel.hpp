/*
 * StreamingMod1DPseudoKernel.cpp
 *
 *  Created on: Apr 10, 2015
 *      Author: pfandedd
 */

#pragma once

#include <vector>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace datadriven {

    template<typename T>
    class StreamingMod1DPseudoKernel {
      private:
        //l == 1 && i == 1 case requires no level and index to be stored
        //as this case leads to a multiplication with 1, the basis functions can be simply filtered out

        //left and right border require no index to be stored (as it is 1)
        std::vector<size_t> gridIndexLeftBorder;
        std::vector<size_t> dim1DLeftBorder;
        std::vector<T> levelLeftBorder;

        std::vector<size_t> gridIndexRightBorder;
        std::vector<size_t> dim1DRightBorder;
        std::vector<T> indexRightBorder;
        std::vector<T> levelRightBorder;

        //base case
        std::vector<size_t> gridIndexBase;
        std::vector<size_t> dim1DBase;
        std::vector<T> levelBase;
        std::vector<T> indexBase;

        T* kernelDataset;
        size_t datasetSize;
        size_t gridSize;

        size_t dim;

        void create1DStructure(T* level, T* index) {
          size_t curDim = 0;
          //iterate over the grid to classify the 1D basis functions
          size_t gridIndex = 0;

          for (size_t basisFunctionIndex = 0;
               basisFunctionIndex < this->gridSize * this->dim;
               basisFunctionIndex++) {
            T l = level[basisFunctionIndex];
            T i = index[basisFunctionIndex];

            if (l == 2.0) {
              //nothing to do
            } else if (i == 1.0) {
              gridIndexLeftBorder.push_back(gridIndex);
              levelLeftBorder.push_back(l);
              dim1DLeftBorder.push_back(curDim);
            } else if (i == l - 1.0) { //l is actual 2^l
              gridIndexRightBorder.push_back(gridIndex);
              levelRightBorder.push_back(l);
              indexRightBorder.push_back(i);
              dim1DRightBorder.push_back(curDim);
            } else {
              gridIndexBase.push_back(gridIndex);
              levelBase.push_back(l);
              indexBase.push_back(i);
              dim1DBase.push_back(curDim);
            }

            curDim = (curDim + 1) % this->dim;

            if (curDim == 0) {
              gridIndex += 1;
            }
          }

        }
      public:
        StreamingMod1DPseudoKernel(size_t dim, T* level, T* index, size_t gridSize,
                                   T* kernelDataset, size_t datasetSize) {
          this->dim = dim;
          this->kernelDataset = kernelDataset;
          this->datasetSize = datasetSize;
          this->gridSize = gridSize;

          create1DStructure(level, index);
        }

        void mult(T* alphaArray, T* result) {
          //      size_t blockSize = 1;
          size_t gridBlockSize = 64;


          //TODO add OMP here
          #pragma omp parallel for

          for (size_t dataIndex = 0; dataIndex < this->datasetSize; dataIndex +=
                 1) {
            //T *dataTuple = this->kernelDataset + (dataIndex * this->dim);

            std::vector<T> dataTuple(dim); //necessary due to transposed dataset

            for (size_t d = 0; d < this->dim; d++) {
              dataTuple[d] = this->kernelDataset[dataIndex
                                                 + (d * this->datasetSize)];
              //        std::cout << "data d: " << d << " -> " << dataTuple[d] << std::endl;
            }

            std::vector<T> eval1DArray(gridBlockSize); //one temporary for each grid point

            T resultComponent = 0.0;

            size_t baseIndex = 0;
            size_t leftIndex = 0;
            size_t rightIndex = 0;

            size_t endBase = this->gridIndexBase.size();
            size_t endLeft = this->gridIndexLeftBorder.size();
            size_t endRight = this->gridIndexRightBorder.size();

            for (size_t gridBlockIndex = 0; gridBlockIndex < this->gridSize;
                 gridBlockIndex += gridBlockSize) {

              //initiliaze with alpha
              for (size_t i = 0; i < gridBlockSize; i++) {
                eval1DArray[i] = alphaArray[gridBlockIndex + i];
              }

              while (baseIndex < this->gridIndexBase.size()) {
                size_t gridPointIndex = this->gridIndexBase[baseIndex];

                if (gridPointIndex >= gridBlockIndex + gridBlockSize) {
                  break;
                }

                T l = this->levelBase[baseIndex];
                T i = this->indexBase[baseIndex];
                size_t curDim = this->dim1DBase[baseIndex];
                T eval1D = std::max(1.0 - fabs(l * dataTuple[curDim] - i),
                                    0.0);
                eval1DArray[gridPointIndex % gridBlockSize] *= eval1D;
                baseIndex += 1;
              }

              //do case left border
              //        for (size_t cur = 0; cur < this->gridIndexLeftBorder.size();
              //            cur++) {
              while (leftIndex < this->gridIndexLeftBorder.size()) {
                size_t gridPointIndex = this->gridIndexLeftBorder[leftIndex];

                if (gridPointIndex >= gridBlockIndex + gridBlockSize) {
                  break;
                }

                T l = this->levelLeftBorder[leftIndex];
                size_t curDim = this->dim1DLeftBorder[leftIndex];
                T eval1D = std::max(2.0 - (l * dataTuple[curDim]), 0.0);
                eval1DArray[gridPointIndex % gridBlockSize] *= eval1D;
                leftIndex += 1;
              }

              //do case right border
              //        for (size_t cur = 0; cur < this->gridIndexRightBorder.size();
              //            cur++) {
              while (rightIndex < this->gridIndexRightBorder.size()) {
                size_t gridPointIndex =
                  this->gridIndexRightBorder[rightIndex];

                if (gridPointIndex >= gridBlockIndex + gridBlockSize) {
                  break;
                }

                T l = this->levelRightBorder[rightIndex];
                T i = this->indexRightBorder[rightIndex];
                size_t curDim = this->dim1DRightBorder[rightIndex];
                T eval1D = std::max(l * dataTuple[curDim] - i + 1.0, 0.0);
                eval1DArray[gridPointIndex % gridBlockSize] *= eval1D;
                rightIndex += 1;
              }

              //do reduce and store intermediate result

              for (size_t i = 0; i < gridBlockSize; i++) {
                resultComponent += eval1DArray[i];
              }
            }

            result[dataIndex] += resultComponent;
          }
        }

    }
    ;

  }
}
