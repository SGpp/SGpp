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
#include <numeric>

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
            // return value
            double score = 0.0;

            // The largest and the second largest PDF of this grip point
            double max = 0.0;
            double second_max = 0.0;

            // Difference between the largest and the second largest PDF of this grip point
            double gridClassDiffs = 0.0;

            // Evaluations of all grids at (param) seq
            std::vector<double> gridEvals;

            // The most and second dominant class number
            size_t max_class = 0;
            size_t second_class = 1;

            // Get the evaluations of seq at all GridSave
            base::DataVector p(storage.getDimension());

            storage.getPoint(seq).getStandardCoordinates(p);

            if (pre_compute) {
                for (size_t i = 0; i < grids.size(); i++) {
                    std::string key = p.toString();
                    gridEvals.push_back(pre_comp_evals.at(i).at(key));
                }
            }
            else {
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
                max_class = 1;
                second_class = 0;
            }

            for (size_t i = 2; i < gridEvals.size(); i++){
                // use >= n not just > as max and second_max can hav same value. Ex:{1,2,3,3}
                if(gridEvals.at(i) >= max){
                    second_max = max;
                    second_class = max_class;
                    max = gridEvals.at(i);
                    max_class = i;
                }
                else if(gridEvals.at(i) > second_max){
                    second_max = gridEvals.at(i);
                    second_class = i;
                }
            }

            gridClassDiffs = max - second_max;

            //Compare to the neighbors

            base::HashGridPoint& gp = storage.getPoint(seq);

            // Find the geometric neighbors
            base::HashGridIterator iter(storage);
            std::vector<size_t> neighSeq;

            for (size_t d = 0; d < storage.getDimension(); d++) {
                // Left neighbor
                if (hasChild(gp, d, true)) {
                    base::HashGridPoint child = base::HashGridPoint(gp);
                    base::HashGridPoint down = base::HashGridPoint(gp);
                    getChild(gp, d, true, child);
                    goDown(child, down, d, false);
                    neighSeq.push_back(storage.getSequenceNumber(down));
                } else {
                    // Check if left-most grid-point on level, which has no left neigh
                    if (gp.getIndex(d) > 1) {
                        base::HashGridPoint gp_c = base::HashGridPoint(gp);
                        base::HashGridPoint up = base::HashGridPoint(gp);
                        goUp(gp_c, up, d, true);
                        neighSeq.push_back(storage.getSequenceNumber(up));
                    }
                }
                // Right neighbor
                if (hasChild(gp, d, false)) {
                    base::HashGridPoint child = base::HashGridPoint(gp);
                    base::HashGridPoint down = base::HashGridPoint(gp);
                    getChild(gp, d, false, child);
                    goDown(child, down, d, true);
                    neighSeq.push_back(storage.getSequenceNumber(down));
                } else {
                    // Check if right-most grid point on level, which has no right neigh
                    if (gp.getIndex(d) < pow(2.0, gp.getLevel(d)) - 1) {
                        base::HashGridPoint gp_c = base::HashGridPoint(gp);
                        base::HashGridPoint up = base::HashGridPoint(gp);
                        goUp(gp_c, up, d, false);
                        neighSeq.push_back(storage.getSequenceNumber(up));
                    }
                }
            }
            // End of finding the geometric neighbors

            // The diffenrence value for each neighbor
            std::vector<double> neighborDiffs;
            // The distance between each neighbor and the current grid point
            std::vector<double> neighborDists;

            // Iterate over all neighbors
            for (size_t i = 0; i < neighSeq.size(); i++) {
                // For each neighbor
                base::HashGridPoint& neighbor = storage.getPoint(neighSeq.at(i));

                std::unique_ptr<base::OperationEval>
                opEval1(op_factory::createOperationEval(*grids.at(max_class)));
                double neighbor_max = opEval1->eval(*alphas.at(max_class), neighSeq.at(i));

                std::unique_ptr<base::OperationEval>
                opEval2(op_factory::createOperationEval(*grids.at(second_class)));
                double neighbor_sec = opEval2->eval(*alphas.at(second_class), neighSeq.at(i));

                neighborDiffs.push_back(neighbor_max-neighbor_sec);
                neighborDists.push_back(getDistance(gp,neighbor));

            }

            double total = 0.0;
            // Calculate the score
            for (unsigned i=0; i<neighborDists.size(); i++) {
                total = total + 1/neighborDists.at(i);
            }

            for (unsigned i=0; i<neighborDiffs.size(); i++) {
                score = score + neighborDiffs.at(i)/(neighborDists.at(i)*total);
            }

            score = score + gridClassDiffs;

            return score;
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

        // Used to get geom. neighbor of non-leaf grid points
        void GridPointBasedCoarseningFunctor::goDown(base::HashGridPoint& gp,
                                                   base::HashGridPoint& down,
                                                   size_t d,
                                                   bool left) const {
            // Child in direction exists? If not stop.
            if (hasChild(gp, d, left)) {
                getChild(gp, d, left, gp);
                goDown(gp, down, d, left);
            } else {
                unsigned int l = gp.getLevel(d);
                unsigned int i = gp.getIndex(d);
                down.set(d, l, i);
            }
        }

        // Use to get geom. neightbor of leaf grid point
        void GridPointBasedCoarseningFunctor::goUp(base::HashGridPoint& gp,
                                                 base::HashGridPoint& up,
                                                 size_t d,
                                                 bool left) const {
            if (isLeftChild(gp, d) != left) {
                getParent(gp, d, up);
            } else {
                getParent(gp, d, gp);
                goUp(gp, up, d, left);
            }
        }

        bool GridPointBasedCoarseningFunctor::hasChild(const base::HashGridPoint& gp,
                                                     size_t d,
                                                     bool left) const {
            base::HashGridIterator iter(grids.at(current_grid_index)->getStorage());
            iter.set(gp);
            if (left) {
                return iter.hintLeft(d);
            } else {
                return iter.hintRight(d);
            }
        }

        bool GridPointBasedCoarseningFunctor::isLeftChild(const base::HashGridPoint& gp,
                                                        size_t d) const {
            size_t i = gp.getIndex(d);
            return ((i + 1) / 2) % 2 == 1;
        }

        void GridPointBasedCoarseningFunctor::getChild(const base::HashGridPoint& gp,
                                                     size_t d, bool left,
                                                     base::HashGridPoint& child) const {
            unsigned int i = gp.getIndex(d);
            unsigned int l = gp.getLevel(d);
            unsigned int new_l = l + 1;
            unsigned int new_i = i * 2;
            if (left) {
                new_i -= 1;
            } else {
                new_i += 1;
            }
            child.set(d, new_l, new_i);
        }

        void GridPointBasedCoarseningFunctor::getParent(const base::HashGridPoint& gp,
                                                      size_t d,
                                                      base::HashGridPoint& par) const {
            unsigned int i = gp.getIndex(d);
            unsigned int l = gp.getLevel(d);
            unsigned int new_l = l - 1;
            unsigned int new_i = (i + 1) / 2;
            if (new_i % 2 == 0) {
                new_i -= 1;
            }
            par.set(d, new_l, new_i);
        }
        double GridPointBasedCoarseningFunctor::getDistance(base::HashGridPoint& gp1,
                                                            base::HashGridPoint& gp2)
        const {
            base::DataVector dist(gp1->getDimension());

            for (size_t d = 0; d < gp1->getDimension(); d++) {
                double coord1 = gp1.getStandardCoordinate(d);
                double coord2 = gp2.getStandardCoordinate(d);
                double diff_sq = coord1-coord2;
                dist.set(d, diff_sq);
            }

            double distance = dist.l2Norm();
            return distance;
        }

        bool GridPointBasedCoarseningFunctor::isWithinSupport(base::HashGridPoint& gp,
                                                         base::DataVector& point)
        const {
            for (size_t d = 0; d < point.getSize(); d++) {
                double coord = gp.getStandardCoordinate(d);
                size_t level = gp.getLevel(d);
                double step = 1.0 / pow(2.0, static_cast<double>(level));
                double lower = coord - step;
                double upper = coord + step;
                if (point.get(d) < lower || point.get(d) > upper) {
                    return false;
                }
            }
            return true;
        }

    }  // namespace datadriven
}  // namespace sgpp
