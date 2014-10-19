#include <cstdlib>
#include <iostream>
#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"

int main(int argc, char **args)
{

    using namespace sg::base;
    using namespace std;

    size_t DIM = 1;
    size_t LEVEL = 2;

    Grid* grid = Grid::createLinearGrid(DIM);
    GridGenerator* grid_gen = grid->createGridGenerator();
    grid_gen->regular(LEVEL);

    double xs[] = {0.3};

    DataMatrix trainData = DataMatrix(xs, 1, 1);

    DataVector classes = DataVector(1);
    classes.set(0, 1);

    DataVector errors = DataVector(1);
    errors.set(0, 3);

    HashRefinement hash_ref = HashRefinement();
    OnlinePredictiveRefinementDimension strategy = OnlinePredictiveRefinementDimension(&hash_ref);
    strategy.setTrainDataset(&trainData);
    strategy.setClasses(&classes);
    strategy.setErrors(&errors);

    OnlinePredictiveRefinementDimension::refinement_map result = OnlinePredictiveRefinementDimension::refinement_map();
    strategy.collectRefinablePoints(grid->getStorage(), 10, &result);

}
