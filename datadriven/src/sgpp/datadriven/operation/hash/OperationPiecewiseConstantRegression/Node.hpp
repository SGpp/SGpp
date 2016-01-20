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
namespace PiecewiseConstantRegression {

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
            base::DataVector &values, bool verbose = false);

    std::vector<size_t> getSupportIndizes(std::vector<float_t> &x, std::vector<float_t> &h,
            std::vector<size_t> &parentSupport);

    float_t getAverage(std::vector<size_t> &supportIndizes, std::vector<float_t> &x, std::vector<float_t> &h);

    float_t getMSE(std::vector<size_t> &supportIndizes, std::vector<float_t> &x, std::vector<float_t> &h,
            float_t supportValue);

    std::unique_ptr<Node> hierarchizeChild(std::vector<float_t> &x, std::vector<float_t> &h,
            std::vector<size_t> &parentSupport, float_t supportValue, float_t targetMSE, size_t targetMaxLevel,
            size_t nextDim, size_t levelLimit);

    void hierarchize(float_t targetMSE, size_t targetMaxLevel, float_t parentValue = 0.0, size_t refineDim = 0,
            size_t levelLimit = 0);

    std::vector<float_t> getLeftChildX(size_t d);

    std::vector<float_t> getRightChildX(size_t d);

    std::vector<float_t> getChildH(size_t d);

    float_t evaluate(std::vector<float_t> &point);

    float_t integrate(SGPP::base::GridIndex &gridPoint, size_t &integratedNodes, size_t levelLimit = 0);

    uint64_t getChildCount();

    uint64_t getHierarchizationMaxLevel();

    float_t getMSE() {
        float_t mse = 0.0;
        for (size_t i = 0; i < dataset.getNrows(); i++) {
            std::vector<float_t> point;
            dataset.getRow(i, point);
            float_t eval = this->evaluate(point);
            mse += (eval - values[i]) * (eval - values[i]);
        }
        return mse;
    }

};

}
}
}
