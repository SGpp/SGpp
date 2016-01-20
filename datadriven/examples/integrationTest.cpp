/*
 * regressionByInterpolation.cpp
 *
 *  Created on: Dec 17, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <random>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/PiecewiseConstantSmoothedRegressionSystemMatrix.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/OperationPiecewiseConstantRegression.hpp>

int main(int argc, char **argv) {

    //    std::string fileName("parabel2d.arff");
    //    std::string fileName("chess_02D_tr.dat.arff");
    //    std::string fileName("friedman_4d.arff");
    //    std::string fileName("chess_05D_3fields_tr.dat.arff");
    //    std::string fileName("chess_3d.arff");
    std::string fileName("parabola_simple_4d.arff");

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset arffDataset = arffTools.readARFF(fileName);

    SGPP::base::DataMatrix &dataset = arffDataset.getTrainingData();
    SGPP::base::DataVector &values = arffDataset.getClasses();

    int maxLevel = 8;
    float_t lambda = 0.0;

    SGPP::datadriven::OperationPiecewiseConstantRegression piecewiseRegressor(dataset, values);

    std::unique_ptr<SGPP::datadriven::PiecewiseConstantRegression::Node> node = piecewiseRegressor.hierarchize(0.001, 20);

    auto grid = std::shared_ptr<SGPP::base::Grid>(SGPP::base::Grid::createLinearGrid(arffDataset.getDimension()));

    auto generator = std::shared_ptr<SGPP::base::GridGenerator>(grid->createGridGenerator());
    generator->regular(maxLevel);

    SGPP::base::OperationMatrix *C = SGPP::op_factory::createOperationLaplace(*grid);
    SGPP::datadriven::PiecewiseConstantSmoothedRegressionSystemMatrix SMatrix(*node, *grid, *C, lambda);

    SGPP::base::DataVector rhs(grid->getStorage()->size());
    for (size_t i = 0; i < 5; i++) {
        SMatrix.generateb(rhs);
    }

    std::cout << "all done!" << std::endl;
    return 0;
}

