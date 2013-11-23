/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef GREEDYPQREFINEMENT_HPP_
#define GREEDYPQREFINEMENT_HPP_

#include "RefinementDecorator.hpp"
#include "../../Grid.hpp"
#include "dataStructures/SortedGridObjectContainer.hpp"
#include <map>

using namespace std;


namespace sg {
namespace base {

typedef map<const char*,size_t> BasisTypes;

class GreedyPQRefinement: public sg::base::RefinementDecorator {
public:
	GreedyPQRefinement(AbstractRefinement* refinement,
			DataMatrix* dataSet,
			DataVector* errors,
			Grid* grid);
	virtual ~GreedyPQRefinement();


	/**
	 * Refines a grid according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	void free_refine(GridStorage* storage, RefinementFunctor* functor);

	/**
	 * Refines a grid by adding additional Subspaces according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	void freeRefineSubspace(GridStorage* storage, RefinementFunctor* functor);

	//Getters and Setters for private Variables.
	const DataMatrix*& getDataSet() const;
	void setDataSet(const DataMatrix*& dataSet);
	const SortedGridObjectContainer& getGridObjectsSortedByError() const;

protected:


	void collectAllRefinablePointsSorted(GridStorage* storage,
										 RefinementFunctor* functor,
										 SortedGridObjectContainer* sortedRefinablePoints);


private:
	SortedGridObjectContainer gridObjectsSortedByError;
	DataMatrix* dataSet;
	BasisTypes basisTypes;
	size_t gridType;
	DataVector* errors;


	void calculateErrorIndicatorForObject(GridStorage* storage, RefinementFunctor* functor, GridPointErrorContainer* refinableObject);
	double basisFunctionEvalHelper(level_t level, index_t index, double value);
	void addErrorObject();

};

} /* namespace base */
} /* namespace sg */
#endif /* GREEDYPQREFINEMENT_HPP_ */

