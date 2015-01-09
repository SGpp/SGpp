#include "../../OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp"
#include "datadriven/operation/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp"

using namespace sg::base;

namespace sg {
namespace datadriven {

OperationMultipleEvalSubspaceCombined::OperationMultipleEvalSubspaceCombined(Grid &grid, DataMatrix &dataset) :
		AbstractOperationMultipleEvalSubspace(grid, dataset) {
	this->paddedDataset = this->padDataset(dataset);
	this->storage = grid.getStorage();
	//this->dataset = dataset;
	this->dim = dataset.getNcols();

	//this->subspaceSize = (2 * this->dim) + 8;
	this->maxGridPointsOnLevel = 0;

#ifdef X86COMBINED_WRITE_STATS
	string prefix("results/data/stats_");
	string fileName(X86COMBINED_WRITE_STATS);
	this->statsFile.open(prefix + fileName, ios::out);;

	this->statsFile << "# name: " << X86COMBINED_WRITE_STATS_NAME << endl;
	this->statsFile << "refinementStep & ";
	this->statsFile << "nonVirtualGridPoints & ";
	this->statsFile << "totalRegularGridPoints & ";
	this->statsFile << "actualGridPoints & ";
	this->statsFile << "largestArraySubspace & ";
	this->statsFile << "largestListSubspace & ";
	this->statsFile << "numberOfListSubspaces & ";
	this->statsFile << "subspaceCount & ";
	this->statsFile << "avrPointsPerSubspace & ";
	this->statsFile << "memoryEstimate & ";
	this->statsFile << "memoryEfficiency";
	this->statsFile << endl;

#endif
}

OperationMultipleEvalSubspaceCombined::~OperationMultipleEvalSubspaceCombined() {
#ifdef X86COMBINED_WRITE_STATS
	this->statsFile.close();
#endif
}

void OperationMultipleEvalSubspaceCombined::prepare() {
	this->allLevelsIndexMap.clear();
	this->allSubspaceNodes.clear();
	this->prepareSubspaceIterator();
}

void OperationMultipleEvalSubspaceCombined::setCoefficients(DataVector &surplusVector) {
	DataVector level(dim);
	DataVector maxIndex(dim);
	DataVector index(dim);

	//TODO: use appropriate types here
	unsigned int curLevel;
	unsigned int curIndex;
	for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
		sg::base::GridIndex *point = this->storage->get(gridIndex);
		for (size_t d = 0; d < this->dim; d++) {
			point->get(d, curLevel, curIndex);
			level.set(d, curLevel);
			index.set(d, curIndex);
			maxIndex.set(d, 1 << curLevel);
		}

		this->setSurplus(level, maxIndex, index, surplusVector.get(gridIndex));
	}
}

//writes a result vector in the order of the points in the grid storage
void OperationMultipleEvalSubspaceCombined::unflatten(DataVector &result) {
	DataVector level(dim);
	DataVector maxIndex(dim);
	DataVector index(dim);

	unsigned int curLevel;
	unsigned int curIndex;
	for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
		sg::base::GridIndex *point = this->storage->get(gridIndex);
		for (size_t d = 0; d < this->dim; d++) {
			point->get(d, curLevel, curIndex);
			level.set(d, curLevel);
			index.set(d, curIndex);
			maxIndex.set(d, 1 << curLevel);
		}
		double surplus;
		bool isVirtual;
		this->getSurplus(level, maxIndex, index, surplus, isVirtual);

		result.set(gridIndex, surplus);
	}
}

void OperationMultipleEvalSubspaceCombined::setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index,
		double value) {
	uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
	uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
	uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
	X86CombinedSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
	subspace.setSurplus(indexFlat, value);
}

void OperationMultipleEvalSubspaceCombined::getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index,
		double &value, bool &isVirtual) {
	uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
	uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
	uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
	X86CombinedSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
	value = subspace.getSurplus(indexFlat);
	if (std::isnan(value)) {
		isVirtual = true;
	} else {
		isVirtual = false;
	}
}

uint32_t OperationMultipleEvalSubspaceCombined::flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index) {
	uint32_t indexFlat = static_cast<uint32_t>(index.get(0)); //TODO ugly conversion
	indexFlat >>= 1;
	for (size_t i = 1; i < dim; i++) {
		int actualDirectionGridPoints = static_cast<int>(maxIndices.get(i)); //TODO ugly conversion
		actualDirectionGridPoints >>= 1;
		indexFlat *= actualDirectionGridPoints;
		uint32_t actualIndex = static_cast<uint32_t>(index.get(i)); //TODO ugly conversion
		actualIndex >>= 1; //divide index by 2, skip even indices
		indexFlat += actualIndex;
	}
	return indexFlat;
}

uint32_t OperationMultipleEvalSubspaceCombined::flattenLevel(size_t dim, size_t maxLevel, DataVector &level) {
	uint32_t levelFlat = 0;
	levelFlat += static_cast<uint32_t>(level.get(dim - 1)); //TODO ugly conversion
	// loop terminates at -1
	for (int i = static_cast<int>(dim - 2); i >= 0; i--) {
		levelFlat *= static_cast<uint32_t>(maxLevel);
		levelFlat += static_cast<int>(level.get(i));
	}
	return levelFlat;
}

DataMatrix *OperationMultipleEvalSubspaceCombined::padDataset(sg::base::DataMatrix &dataset) {
	size_t chunkSize = X86COMBINED_PARALLEL_DATA_POINTS;

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = dataset.getNrows() % chunkSize;
	size_t loopCount = chunkSize - remainder;

	if (loopCount == chunkSize) {
		return &dataset;
	}

	sg::base::DataVector lastRow(dataset.getNcols());
	size_t oldSize = dataset.getNrows();
	dataset.getRow(dataset.getNrows() - 1, lastRow);

	DataMatrix *paddedDataset = new DataMatrix(dataset);
	//pad to make dataset % X86COMBINED_PARALLEL_DATA_POINTS == 0
	paddedDataset->resize(dataset.getNrows() + loopCount);

	for (size_t i = 0; i < loopCount; i++) {
		paddedDataset->setRow(oldSize + i, lastRow);
	}

	//additional padding for subspace skipping
	//if validIndices contain X86COMBINED_PARALLEL_DATA_POINTS - 1 it is possible for a vector iteration to contain
	//indices larger than size(dataset) (even though the dataset is divided by X86COMBINED_PARALLEL_DATA_POINTS)
	//add X86COMBINED_VEC_PADDING dummy data points to avoid that problem
	//add X86COMBINED_VEC_PADDING * 2 to also enable the calculateIndexCombined2() method
	paddedDataset->addSize(X86COMBINED_VEC_PADDING * 2);
	for (size_t i = paddedDataset->getNrows(); i < paddedDataset->getNrows() + paddedDataset->getUnused(); i++) {
		for (size_t j = 0; j < paddedDataset->getNcols(); j++) {
			paddedDataset->set(i, j, 0.0);
		}
	}

	return paddedDataset;

	//TODO initialize padding area - hacky -> remove this? -> should not be required as long as vector size divides chunksize
//	this->paddedDataset.addSize(X86COMBINED_VEC_PADDING * 2);
//	for (size_t i = this->dataset.getNrows(); i < this->dataset.getNrows() + this->dataset.getUnused(); i++) {
//		for (size_t j = 0; j < this->dataset.getNcols(); j++) {
//			this->dataset.set(i, j, 0.0);
//		}
//	}
}

size_t OperationMultipleEvalSubspaceCombined::getPaddedDatasetSize() {
	return this->paddedDataset->getNrows();
}

size_t OperationMultipleEvalSubspaceCombined::getAlignment() {
	return X86COMBINED_PARALLEL_DATA_POINTS;
}

std::string OperationMultipleEvalSubspaceCombined::getImplementationName() {
	return "COMBINED";
}

}
}
