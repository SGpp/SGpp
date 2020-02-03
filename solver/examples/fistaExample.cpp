// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_fistaExample_cpp FISTA Solver
 * This example demonstrates the FISTA solver for a toy dataset using
 * using the elastic net regularization method with various regularization
 * penalties.
 */

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/solver/sle/fista/ZeroFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/ElasticNetFunction.hpp>

#include <cmath>
#include <vector>
#include <random>

int main(int argc, char** argv) {
    using sgpp::base::DataVector;
    using sgpp::base::DataMatrix;
    /**
     * Create a two-dimensional grid with level 6.
     */
    auto grid = sgpp::base::Grid::createModLinearGrid(2);
    grid->getGenerator().regular(6);

    /**
     * Then create a two dimensional dataset.
     */
    const auto num_examples = 3000;
    auto dataset = DataMatrix(num_examples, 2);
    auto y = DataVector(num_examples);
    auto gen = std::mt19937_64(42);
    auto dist = std::normal_distribution<double>(0, 0.1);
    for (auto i = 0; i < num_examples; ++i) {
        const double x1 = std::abs(std::sin(i));
        const double x2 = std::abs(std::sin(i*x1));
        const double comb = std::sinh(x1) * (x1 + x2) + dist(gen);
        const auto row = DataVector(std::vector<double>({x1, x2}));
        dataset.setRow(i, row);
        y.set(i, comb);
    }

    std::cout << "Created grid with size: " << grid->getSize() << std::endl;
    auto op = sgpp::op_factory::createOperationMultipleEval(
                *grid, dataset);
    /**
     * Set up the weights for the grid.
     */
    DataVector weights(grid->getSize());
    double lambda = 1.0;
    const double gamma = 0.98;
    /**
     * Iterate over a few regularization penalties.
     */
    for (int i = 0; i < 5; ++i) {
        // weights.setAll(0.0);
        lambda /= 10;
    /**
    * Create an elastic net function and the corresponding fista solver.
    */
        auto g = sgpp::solver::ElasticNetFunction(lambda, gamma);
        auto solver = sgpp::solver::Fista<decltype(g)>(g);
        std::cout << "Start solving." << std::endl;
    /**
    * Here we solve the penalised linear system.
    */
        solver.solve(*op, weights, y, 500, 1e-4);
        std::cout << "Finished solving." << std::endl;
        auto prediction = DataVector(num_examples);
    /**
    * Finally, we calculate the root-mean-squared error of the prediction
    * and print it.
    */
        op->mult(weights, prediction);
        prediction.sub(y);
        prediction.sqr();
        const double rmse = std::sqrt((1.0/num_examples) * prediction.sum());
        std::cout << "Lambda = " << lambda
                  << " |weights|_2 = " << weights.l2Norm() << std::endl;
        std::cout << "Residual = " << rmse << std::endl;
    }
}
