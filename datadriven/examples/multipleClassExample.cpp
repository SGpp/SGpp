// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>


#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sgpp/datadriven/application/MultipleClassPoint.hpp>
#include <sgpp/datadriven/tools/SimpleMultiClassGenerator.hpp>


/**
 * TODO
 * \page example_classificationRefinementExample_cpp Classification Example
 *
 * This example shows how classification specific refinement strategies
 * are used. To do classification, for each class a PDF is approximated with
 * LearnerSGDE and the class with the highest probability gets assigned for
 * new data points to be classified.
 * The ripley data sets is used, although the small number of training data
 * poitns in combination with only a basic setup does not yield good results
 * for any refinement strategy. This example is merely a tech-example.
 */


int main() {
    std::cout << std::endl << "Start" << std::endl;
    size_t dim = 2;
    int classes = 3;
    int level = 3;

    // TODO what to set/where when?
    std::vector<sgpp::base::DataVector*> alphas;
    size_t numRefinements = 3;
    bool levelPenalize = false;
    bool preCompute = false;
    double thresh = 0.0;

    sgpp::datadriven::SimpleMultiClassGenerator generator =
            sgpp::datadriven::SimpleMultiClassGenerator(dim, classes, level);
    std::vector<sgpp::base::Grid*> grids = generator.getMulitGrid();

    std::vector<sgpp::datadriven::MultipleClassPoint> points;

    std::cout << "Size: " << grids.at(0)->getSize() << std::endl;
    std::cout << "Dim: " << grids.at(0)->getDimension() << std::endl;
    std::cout << std::endl;

    for (size_t i = 0; i < grids.at(0)->getStorage().getSize(); i++) {
        //sgpp::base::GridPoint& gp = grids.at(0)->getStorage().getPoint(i);
        points.push_back(sgpp::datadriven::MultipleClassPoint(classes, i, generator));

        //std::cout << "Point " << i << ": " << points.at(i).getDominateClass();
        //std::cout << " ~ " << gp.getStandardCoordinate(0) << " - " << gp.getStandardCoordinate(1);
        //std::cout << std::endl <<  points.at(i).toString() << std::endl;
      }

    sgpp::datadriven::MultipleClassRefinementFunctor mcrf(grids, alphas, &points,
            grids.at(0)->getStorage(), numRefinements, levelPenalize, preCompute, thresh);
    sgpp::datadriven::MultiGridRefinementFunctor* mulitfun = &mcrf;
    
    double totalScore = mulitfun->getTotalRefinementValue(grids.at(0)->getStorage());
    std::cout << "TotalRefinementValue: " << totalScore << std::endl;

    std::cout << std::endl << std::endl;
    for (size_t i = 0; i < grids.at(0)->getStorage().getSize(); i++) {
        sgpp::base::GridPoint& gp = grids.at(0)->getStorage().getPoint(i);

        std::cout << "Point " << i << ": " << points.at(i).getDominateClass();
        std::cout << " ~ " << gp.getStandardCoordinate(0) << " - " << gp.getStandardCoordinate(1);
        std::cout << std::endl <<  points.at(i).toString() << std::endl;
      }

    std::cout << std::endl << "Done" << std::endl;
    return 0;

}
