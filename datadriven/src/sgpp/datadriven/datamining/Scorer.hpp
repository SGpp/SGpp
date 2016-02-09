/*
 * Scorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SCORER_HPP_
#define SCORER_HPP_


#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/datamining/Metric.hpp>
#include <sgpp/datadriven/datamining/ModelFittingBase.hpp>

#include <memory>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * Superclass for TestError, TrainError, CrossValidation, AIC, etc.
 */
class Scorer {
public:
	Scorer(std::shared_ptr<Metric> metric, std::shared_ptr<ModelFittingBase> fitter):metric(metric), fitter(fitter){};
	virtual ~Scorer();
	virtual double getScore(const Dataset& dataset) = 0;

protected:
	std::shared_ptr<Metric> metric;
	std::shared_ptr<ModelFittingBase> fitter;
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* SCORER_HPP_ */
