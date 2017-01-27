/*
 * SimpleMultiClassGenerator.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: katrin
 */

#include "SimpleMultiClassGenerator.hpp"

#include <sgpp/base/grid/Grid.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
SimpleMultiClassGenerator::SimpleMultiClassGenerator(int dim, int classes, int level)
{
    for (int i = 0 ; i < classes ; i++) {
        sgpp::base::Grid* multigrid(sgpp::base::Grid::createLinearGrid(dim));
        multigrid->getGenerator().regular(level);
        grids.push_back(multigrid);
    }
}

SimpleMultiClassGenerator::~SimpleMultiClassGenerator()
{
    // TODO Auto-generated destructor stub
}
std::vector<base::Grid*> SimpleMultiClassGenerator::getMulitGrid() {
    return grids;
}

double SimpleMultiClassGenerator::getEval(int classId, int seq) {
    base::HashGridPoint& gp = grids.at(classId)->getStorage().getPoint(seq);
    base::DataVector coords(grids.at(classId)->getDimension());
    gp.getStandardCoordinates(coords);
    double evals = 0.0;
    switch (classId) {
        case 0:
            evals = coords.sum() / (double) coords.getSize();
            break;
        case 1:
            evals = 1 - (coords.sum() / (double) coords.getSize());
            break;
        case 2:
            evals = 0.6;
            break;
        case 3:
            for (int i = 0; i < coords.getSize(); i++) {
                evals = evals * coords.get(i);
            }
            break;
        default:
            evals = 0.0;
    }
    return evals;
}

} /* namespace datadriven */
} /* namespace sgpp */