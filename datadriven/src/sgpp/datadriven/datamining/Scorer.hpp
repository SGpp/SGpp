/*
 * Scorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SCORER_HPP_
#define SCORER_HPP_

#include <sgpp/datadriven/datamining/Metric.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * Superclass for TestError, TrainError, CrossValidation, AIC, etc.
 */
class Scorer {
public:
	Scorer(Metric* metric):metric(metric){};
	virtual ~Scorer();

private:
	Metric* metric;
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* SCORER_HPP_ */
