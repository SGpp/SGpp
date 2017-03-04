// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIPLECLASSREFINEMENTFUNCTOR_HPP
#define MULTIPLECLASSREFINEMENTFUNCTOR_HPP

#include "ZeroCrossingRefinementFunctor.hpp"


#include <sgpp/datadriven/application/MultipleClassPoint.hpp>
#include <vector>
#include <tuple>

namespace sgpp {
namespace datadriven {

class MultipleClassRefinementFunctor: public ZeroCrossingRefinementFunctor {
public:
	MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num,
                                bool level_penalize,
                                bool pre_compute,
                                double thresh);
                                
	double operator()(base::GridStorage& storage,
                    size_t seq) const override;
                    
    base::Grid* getCombinedGrid();
                    

private:
	std::vector<sgpp::datadriven::MultipleClassPoint> points;
	base::Grid* multigrid;

	void findCrossings(int leftP, int rightP, int seq, size_t d);
	bool hasChild(const base::HashGridPoint& gp, size_t d, bool left) const;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULTIPLECLASSREFINEMENTFUNCTOR_HPP */
