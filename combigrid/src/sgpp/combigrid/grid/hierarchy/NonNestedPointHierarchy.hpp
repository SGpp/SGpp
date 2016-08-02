/*
 * NonNestedPointHierarchy.hpp
 *
 *  Created on: 18.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_HIERARCHY_NONNESTEDPOINTHIERARCHY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_HIERARCHY_NONNESTEDPOINTHIERARCHY_HPP_

#include <sgpp/globaldef.hpp>
#include "AbstractPointHierarchy.hpp"
#include "../ordering/AbstractPointOrdering.hpp"
#include "../distribution/AbstractPointDistribution.hpp"
#include "../../common/AbstractPermutationIterator.hpp"

#include <vector>

namespace sgpp{
namespace combigrid {

/**
 * PointHierarchy class for point hierarchies that are not fully nested (see also AbstractPointHierarchy).
 */
class NonNestedPointHierarchy: public AbstractPointHierarchy {
public:
	std::vector<size_t> numPointsPerLevel;
	std::vector<std::vector<double>> points;
	std::vector<std::shared_ptr<AbstractPermutationIterator>> permutationIterators;

	std::shared_ptr<AbstractPointDistribution> pointDistribution;
	std::shared_ptr<AbstractPointOrdering> pointOrdering;

public:
	NonNestedPointHierarchy(std::shared_ptr<AbstractPointDistribution> pointDistribution, std::shared_ptr<AbstractPointOrdering> pointOrdering);

	virtual ~NonNestedPointHierarchy();

	virtual double getPoint(size_t level, size_t index);

	virtual std::vector<double> &computePoints(size_t level);
	virtual std::vector<double> getPoints(size_t level, bool sorted);

	virtual size_t getNumPoints(size_t level);

	virtual bool isNested();

	virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_HIERARCHY_NONNESTEDPOINTHIERARCHY_HPP_ */
