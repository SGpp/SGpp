// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

namespace SGPP {
namespace datadriven {
namespace HistogramTree {

class Node {
private:
    std::vector<float_t> x;
    std::vector<float_t> h;

    size_t dim;

    std::vector<size_t> supportIndizes;

    std::unique_ptr<Node> leftChild;
    std::unique_ptr<Node> rightChild;
    size_t childDim;

    float_t surplus;

    base::DataMatrix &dataset;
    base::DataVector &values;

    size_t childCount;

    bool verbose;
public:

    static uint64_t integratedNodes;

    static uint64_t hierarchizeMaxLevel;

    Node(std::vector<float_t> x, std::vector<float_t> h, std::vector<size_t> supportIndizes, base::DataMatrix &dataset,
            base::DataVector &values, bool verbose = false) :
            x(x), h(h), dim(x.size()), supportIndizes(supportIndizes), leftChild(nullptr), rightChild(nullptr), childDim(
                    0), surplus(0.0), dataset(dataset), values(values), childCount(0), verbose(verbose) {
    }

    std::vector<size_t> getSupportIndizes(std::vector<float_t> &x, std::vector<float_t> &h,
            std::vector<size_t> &parentSupport) {
        std::vector<size_t> supportIndizes;

        for (size_t i = 0; i < parentSupport.size(); i++) {
            size_t dataIndex = parentSupport[i];

            base::DataVector point(dim);
            dataset.getRow(dataIndex, point);
            bool onSupport = true;

            for (size_t d = 0; d < dim; d++) {
                if (point[d] < (x[d] - h[d]) || point[d] > (x[d] + h[d])) {
                    onSupport = false;
                    break;
                }
            }

            if (onSupport) {
                supportIndizes.push_back(dataIndex);
            }
        }
        return supportIndizes;
    }

    float_t getAverage(std::vector<size_t> &supportIndizes, std::vector<float_t> &x, std::vector<float_t> &h) {
        float_t sum = 0.0;
        size_t pointOnSupport = 0;

        for (size_t i = 0; i < supportIndizes.size(); i++) {
            size_t dataIndex = supportIndizes[i];
            sum += values[dataIndex];
            pointOnSupport += 1;
        }

        float_t average;
        if (pointOnSupport > 0) {
            average = sum / static_cast<float_t>(pointOnSupport);
        } else {
            average = 0.0;
        }
        return average;
    }

    float_t getMSE(std::vector<size_t> &supportIndizes, std::vector<float_t> &x, std::vector<float_t> &h,
            float_t supportValue) {

        float_t sum = 0.0;
        size_t pointOnSupport = 0;

        for (size_t i = 0; i < supportIndizes.size(); i++) {
            size_t dataIndex = supportIndizes[i];

            float_t diff = supportValue - values[dataIndex];
            sum += diff * diff;
            pointOnSupport += 1;
        }

        float_t mse;
        if (pointOnSupport > 0) {
            mse = sum / static_cast<float_t>(pointOnSupport);
        } else {
            mse = 0.0;
        }
        return mse;
    }

    std::unique_ptr<Node> hierarchizeChild(std::vector<float_t> &x, std::vector<float_t> &h,
            std::vector<size_t> &parentSupport, float_t supportValue, float_t targetMSE, size_t targetMaxLevel,
            size_t nextDim, size_t levelLimit) {
        std::vector<size_t> support = getSupportIndizes(x, h, parentSupport);

        if (support.size() == 0) {
            if (verbose) {
                std::cout << "reached 0-points: " << std::endl;
            }
            return std::unique_ptr<Node>(nullptr);
        }

        float_t mse = getMSE(support, x, h, supportValue);

        if (mse <= targetMSE) {
            if (verbose) {
                std::cout << "reached MSE: " << targetMSE << std::endl;
            }
            return std::unique_ptr<Node>(nullptr);
        }

        std::unique_ptr<Node> childNode = std::make_unique<Node>(x, h, support, dataset, values);
        childNode->hierarchize(targetMSE, targetMaxLevel, supportValue, nextDim, levelLimit);

        return std::move(childNode);
    }

    void hierarchize(float_t targetMSE, size_t targetMaxLevel, float_t parentValue = 0.0, size_t refineDim = 0,
            size_t levelLimit = 0) {

        if (levelLimit > targetMaxLevel) {
            return;
        }

        if (levelLimit > hierarchizeMaxLevel) {
            Node::hierarchizeMaxLevel = levelLimit;
        }

        surplus = getAverage(supportIndizes, x, h) - parentValue;
        float_t supportValue = parentValue + surplus;

//        for (size_t d = refineDim; d < dim; d++) {

        bool refinedAnyChild = false;

        std::vector<float_t> childH = getChildH(refineDim);

        std::vector<float_t> leftChildX = getLeftChildX(refineDim);
        leftChild = hierarchizeChild(leftChildX, childH, supportIndizes, supportValue, targetMSE, targetMaxLevel,
                (refineDim + 1) % dim, levelLimit + 1);
        if (leftChild.operator bool()) {
            refinedAnyChild = true;
            childCount += 1 + leftChild->childCount;
        }

        std::vector<float_t> rightChildX = getRightChildX(refineDim);
        rightChild = hierarchizeChild(rightChildX, childH, supportIndizes, supportValue, targetMSE, targetMaxLevel,
                (refineDim + 1) % dim, levelLimit + 1);
        if (rightChild.operator bool()) {
            refinedAnyChild = true;
            childCount += 1 + rightChild->childCount;
        }

        if (refinedAnyChild) {
            this->childDim = refineDim;
//                break;
        }
//        }

    }

    std::vector<float_t> getLeftChildX(size_t d) {
        std::vector<float_t> childX(x);
        childX[d] -= (h[d] / 2.0);
        return childX;
    }

    std::vector<float_t> getRightChildX(size_t d) {
        std::vector<float_t> childX(x);
        childX[d] += (h[d] / 2.0);
        return childX;
    }

    std::vector<float_t> getChildH(size_t d) {
        std::vector<float_t> childH(h);
        childH[d] = h[d] / 2.0;
        return childH;
    }

    float_t evaluate(std::vector<float_t> &point) {
        float_t sum = 0.0;
//        std::cout << "evaluate at x:";
//        for(size_t d = 0; d < dim; d++) {
//            if (d > 0) {
//                std::cout << ", ";
//            }
//            std::cout << x[d];
//        }
//        std::cout << " h: ";
//        for(size_t d = 0; d < dim; d++) {
//            if (d > 0) {
//                std::cout << ", ";
//            }
//            std::cout << h[d];
//        }
//        std::cout << " surplus: " << surplus;
//        std::cout << std::endl;

        sum += surplus;
        if (point[this->childDim] < x[this->childDim]) {
            if (this->leftChild.operator bool()) {
                sum += this->leftChild->evaluate(point);
            }
        } else {
            if (this->rightChild.operator bool()) {
                sum += this->rightChild->evaluate(point);
            }
        }

        return sum;
    }

//#define LEVEL_TO_PRINT 0
//#define INDEX_TO_PRINT 0

    float_t integrate(SGPP::base::GridIndex &gridPoint, size_t levelLimit = 0) {

        if (levelLimit == 0) {
            Node::integratedNodes = 1;
        } else {
            Node::integratedNodes += 1;
        }

//        if (levelLimit > 2) {
//            return 0.0;
//        }
        SGPP::base::SLinearBase basis;

        float_t sum = 0.0;

//        if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//            std::cout << "x: ";
//            for (size_t d = 0; d < dim; d++) {
//                if (d > 0) {
//                    std::cout << ", ";
//                }
//                std::cout << x[d];
//            }
//
//            std::cout << " h: ";
//            for (size_t d = 0; d < dim; d++) {
//                if (d > 0) {
//                    std::cout << ", ";
//                }
//                std::cout << h[d];
//            }
//            std::cout << std::endl;
//            std::cout << "surplus: " << surplus << std::endl;
//        }

        float_t integral = 1.0;
        for (size_t d = 0; d < dim; d++) {
            //integrate left side of triangle
            float_t integral1D = 0.0;

            float_t gridPointHat = gridPoint.getCoord(d);

            float_t gridPointH = (1.0 / static_cast<float_t>(1 << gridPoint.getLevel(d)));
//                    * static_cast<float_t>(gridPoint.getIndex(d));

//            if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                std::cout << "gridPointH: " << gridPointH << std::endl;
//            }

            float_t leftGridPointHat = gridPointHat - gridPointH;
            float_t rightGridPointHat = gridPointHat + gridPointH;

//            if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                std::cout << "leftGridPointHat: " << leftGridPointHat << std::endl;
//                std::cout << "rightGridPointHat: " << rightGridPointHat << std::endl;
//            }

            float_t leftGridPointConstant = x[d] - h[d];
            float_t rightGridPointConstant = x[d] + h[d];

//            if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                std::cout << "leftGridPointConstant: " << leftGridPointConstant << std::endl;
//                std::cout << "rightGridPointConstant: " << rightGridPointConstant << std::endl;
//            }

            float_t leftSideLeftBorder = std::max(leftGridPointHat, leftGridPointConstant);
            float_t leftSideRightBorder = std::min(gridPointHat, rightGridPointConstant);

            if (leftSideRightBorder > leftSideLeftBorder) {
//                if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                    std::cout << "leftSideLeftBorder: " << leftSideLeftBorder << std::endl;
//                    std::cout << "leftSideRightBorder: " << leftSideRightBorder << std::endl;
//                }

                float_t leftIntegral1D = (leftSideRightBorder - leftSideLeftBorder)
                        * basis.eval(gridPoint.getLevel(d), gridPoint.getIndex(d),
                                (leftSideRightBorder + leftSideLeftBorder) / 2.0);
                integral1D += leftIntegral1D;

//                if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                    std::cout << "left integral: " << leftIntegral1D << std::endl;
//                }
            }

            float_t rightSideLeftBorder = std::max(gridPointHat, leftGridPointConstant);
            float_t rightSideRightBorder = std::min(rightGridPointHat, rightGridPointConstant);

            if (rightSideRightBorder > rightSideLeftBorder) {
//                if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                    std::cout << "rightSideLeftBorder: " << rightSideLeftBorder << std::endl;
//                    std::cout << "rightSideRightBorder: " << rightSideRightBorder << std::endl;
//                }

                float_t rightIntegral1D = (rightSideRightBorder - rightSideLeftBorder)
                        * basis.eval(gridPoint.getLevel(d), gridPoint.getIndex(d),
                                (rightSideRightBorder + rightSideLeftBorder) / 2.0);
                integral1D += rightIntegral1D;

//                if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                    std::cout << "right integral: " << rightIntegral1D << std::endl;
//                }
            }

            integral *= integral1D;

//            if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//                std::cout << "---dim---" << std::endl;
//            }
        }

        float_t product = surplus * integral;

//        if (gridPoint.getLevel(0) == LEVEL_TO_PRINT and gridPoint.getIndex(0) == INDEX_TO_PRINT) {
//            std::cout << "----------------------" << std::endl;
//        }
        if (integral > 0.0) {
            if (leftChild.operator bool()) {
                sum += leftChild->integrate(gridPoint, levelLimit + 1);
            }
            if (rightChild.operator bool()) {
                sum += rightChild->integrate(gridPoint, levelLimit + 1);
            }
        }
        return sum + product;
    }

    uint64_t getChildCount() {
        return childCount;
    }

    uint64_t getHierarchizationMaxLevel() {
        return hierarchizeMaxLevel;
    }

};

}
}
}
