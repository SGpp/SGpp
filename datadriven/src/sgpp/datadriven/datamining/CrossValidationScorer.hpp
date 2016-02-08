/*
 * SimpleSplittingScorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef CROSSVALIDATIONSCORERSCORER_HPP_
#define CROSSVALIDATIONSCORERSCORER_HPP_

#include <vector>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include "Scorer.hpp"

namespace SGPP {
namespace datadriven {

class CrossValidationScorer: public Scorer {
public:
	CrossValidationScorer(Dataset& dataset, Metric& metric, size_t kFold);
	virtual ~CrossValidationScorer();

	void splitset(base::DataMatrix& dataset,
	                            base::DataVector& datasetValues, size_t kFold,
	                            std::vector<base::DataMatrix>& trainingSets,
	                            std::vector<base::DataVector>& trainingSetValues,
	                            std::vector<base::DataMatrix>& testSets,
	                            std::vector<base::DataVector>& testSetValues, bool verbose);

private:
	Dataset trainDataset;
	Dataset testDataset;
	size_t kFold;

};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* CROSSVALIDATIONSCORERSCORER_HPP_ */
