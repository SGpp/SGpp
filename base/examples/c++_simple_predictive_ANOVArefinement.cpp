// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <math.h>
// All SG++ headers
//#include <sgpp_base.hpp>

// Or, better!, include only those that are required
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveANOVARefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

using namespace std;
using namespace SGPP::base;

// function to interpolate
SGPP::float_t f(SGPP::float_t x0, SGPP::float_t x1) {
    return sin(x0 * 10) + x1;
}

DataVector& calculateError(const DataMatrix& dataSet, Grid& grid, const DataVector& alpha, DataVector& error) {
    cout << "calculating error" << endl;

    //traverse dataSet
    DataVector vec(2);
    OperationEval* opEval = SGPP::op_factory::createOperationEval(grid);
    for (uint i = 0; i < dataSet.getNrows(); i++) {
        dataSet.getRow(i, vec);
        error[i] = pow(f(dataSet.get(i, 0), dataSet.get(i, 1)) - opEval->eval(alpha, vec), 2);
    }
    return error;
}

int main() {
    // create a two-dimensional piecewise bilinear grid
    size_t dim = 2;
    Grid* grid = Grid::createModLinearGrid(dim);
    GridStorage* hashGridStorage = grid->getStorage();
    cout << "dimensionality:         " << hashGridStorage->dim() << endl;

    // create regular grid, level 3
    size_t level = 3;
    GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(level);
    cout << "number of grid points:  " << hashGridStorage->size() << endl;

    // create coefficient vector
    DataVector alpha(hashGridStorage->size());
    alpha.setAll(0.0);
    cout << "length of alpha vector: " << alpha.getSize() << endl;

    int rows = 100;
    int cols = 100;

    DataMatrix dataSet(rows * cols, dim);
    DataVector vals(rows * cols);

    // Create a "List" of points where the error should be calculated.
    // This represents a regular 2d grid with a step size of 1 / rows and 1 / cols.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            //xcoord
            dataSet.set(i * cols + j, 0, i * 1.0 / rows);
            //ycoord
            dataSet.set(i * cols + j, 1, j * 1.0 / cols);
            vals[i * cols + j] = f(i * 1.0 / rows, j * 1.0 / cols);
        }
    }

    //refine adaptively 20 times
    for (int step = 0; step < 20; step++) {
        // set function values in alpha
        DataVector gridPointCoordinates(dim);

        for (size_t i = 0; i < hashGridStorage->size(); i++) {
            hashGridStorage->get(i)->getCoords(gridPointCoordinates);
            alpha[i] = f(gridPointCoordinates[0], gridPointCoordinates[1]);
        }

        // hierarchize
        SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);

        // calculate squared offset
        DataVector errorVector(dataSet.getNrows());
        calculateError(dataSet, *grid, alpha, errorVector);

        // refinement  stuff
        HashRefinement refinement;
        PredictiveANOVARefinement decorator(&refinement);

        // refine a single grid point each time
        cout << "Error over all = "  << errorVector.sum() << endl;
        PredictiveRefinementIndicator indicator(grid, &dataSet, &errorVector, 1);
        decorator.free_refine(hashGridStorage, &indicator);

        cout << "Refinement step " << step + 1 << ", new grid size: " << hashGridStorage->size() << endl;

        // plot grid

        // extend alpha vector (new entries uninitialized)
        alpha.resize(hashGridStorage->size());
    }

    delete grid;
}
