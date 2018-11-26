// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedCoarseningFunctor.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>
#include <map>
#include <string>


namespace sgpp {
    namespace datadriven {

        GridPointBasedCoarseningFunctor::
        GridPointBasedCoarseningFunctor(std::vector<base::Grid*> grids,
                                        std::vector<base::DataVector*> alphas,
                                        size_t r_num,
                                        bool pre_compute,
                                        double thresh) :
                grids(grids), alphas(alphas), current_grid_index(0),
                coarsenings_num(r_num), threshold(thresh),
                pre_compute(pre_compute),
                pre_comp_evals() {
            // Add a map for each grid
            for (size_t i = 0; i < grids.size(); i++) {
                pre_comp_evals.push_back(std::map<std::string, double>());
            }
        }

        double
        GridPointBasedCoarseningFunctor::operator()(base::GridStorage& storage,
                                                    size_t seq) const {

            // The largest and the second largest PDF of this grip point
            double max = 0.0;
            double second_max = 0.0;

            // Difference between the largest and the second largest PDF of this grip point
            double gridClassDiffs = 0.0;

            // Evaluations of all grids at (param) seq
            std::vector<double> gridEvals;



            // Get the evaluations of seq at all GridSave
            base::DataVector p(storage.getDimension());
            storage.getPoint(seq).getStandardCoordinates(p);
            if (pre_compute) {
                for (size_t i = 0; i < grids.size(); i++) {
                    std::string key = p.toString();
                    gridEvals.push_back(pre_comp_evals.at(i).at(key));
                }
            } else {
                for (size_t i = 0; i < grids.size(); i++) {
                    std::unique_ptr<base::OperationEval>
                    opEval(op_factory::createOperationEval(*grids.at(i)));
                    gridEvals.push_back(opEval->eval(*alphas.at(i), p));
                }
            }

            // Find the largest and the second largest PDF of this grip point
            if(gridEvals.at(0) > gridEvals.at(1)) {
                second_max = gridEvals.at(1);
                max = gridEvals.at(0);
            } else {
                second_max = gridEvals.at(0);
                max = gridEvals.at(1);
            }

            for (int i = 2; i < gridEvals.size(); i++){
                // use >= n not just > as max and second_max can hav same value. Ex:{1,2,3,3}
                if(gridEvals.at(i) >= max){
                    second_max = max;
                    max = gridEvals.at(i);
                }
                else if(gridEvals.at(i) > second_max){
                    second_max = gridEvals.at(i);
                }
            }

            gridClassDiffs = max - second_max;

            return gridClassDiffs;
        }

        void GridPointBasedCoarseningFunctor::preComputeEvaluations() {
            base::DataVector p(grids.at(0)->getDimension());
            std::string key = "";
            double v = 0.0;
            // Evaluated at (!) grid with index i and store in the i-th map
            // Here grid points are not only evaluated at their own grid but at
            // all grids
            for (size_t i = 0; i < grids.size(); i++) {
                std::unique_ptr<base::OperationEval>
                opEval(op_factory::createOperationEval(*grids.at(i)));
                pre_comp_evals.at(i).clear();
                // Iterate over all possible grid point coordinates of all grids
                for (size_t j = 0; j < grids.size(); j++) {
                    for (size_t k = 0; k < grids.at(j)->getSize(); k++) {
                        // The coordinates
                        grids.at(j)->getStorage().getPoint(k).getStandardCoordinates(p);
                        // The hash key
                        key = p.toString();
                        // Does the key already exist?
                        if (pre_comp_evals.at(i).count(key) == 0) {
                            v = opEval->eval(*alphas.at(i), p);
                            pre_comp_evals.at(i).insert(std::pair<std::string, double>(key,
                                                                                       v));
                        }
                    }
                }
            }
        }

        double GridPointBasedCoarseningFunctor::start() const {
            return 0.0;
        }

        size_t GridPointBasedCoarseningFunctor::getCoarseningsNum() const {
            return this->coarsenings_num;
        }

        double GridPointBasedCoarseningFunctor::getCoarseningThreshold() const {
            return this->threshold;
        }

        void GridPointBasedCoarseningFunctor::setGridIndex(size_t grid_index) {
            this->current_grid_index = grid_index;
        }

        size_t GridPointBasedCoarseningFunctor::getNumGrids() {
            return this->grids.size();
        }

    }  // namespace datadriven
}  // namespace sgpp
