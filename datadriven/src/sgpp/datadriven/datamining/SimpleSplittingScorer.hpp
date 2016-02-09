/*
 * SimpleSplittingScorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SIMPLESPLITTINGSCORER_HPP_
#define SIMPLESPLITTINGSCORER_HPP_

#include <vector>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include "Scorer.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

class SimpleSplittingScorer: public Scorer {
public:
	SimpleSplittingScorer(Dataset& dataset, Metric* metric, double trainPortion, int seed=-1);
	virtual ~SimpleSplittingScorer();

	void splitset(const DataMatrix& dataset,
            const DataVector& datasetValues, double trainPortion,
            DataMatrix& trainingSet,
			 DataVector& trainingSetValues,
			 DataMatrix& testSet,
			 DataVector& testSetValues, bool permute) ;

private:
	Dataset trainDataset;
	Dataset testDataset;
	double trainPortion;
	long int seed;

};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* SIMPLESPLITTINGSCORER_HPP_ */
