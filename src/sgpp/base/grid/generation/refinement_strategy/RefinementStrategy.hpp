/*
 * RefinementStrategy.hpp
 *
 *  Created on: Mar 6, 2012
 *      Author: khakhutv_local
 */

#ifndef REFINEMENTSTRATEGY_HPP_
#define REFINEMENTSTRATEGY_HPP_

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
//#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"


namespace sg
{
namespace base
{

class HashRefinementAbstract;

class RefinementStrategy {
public:
	/*RefinementStrategy();*/
	RefinementStrategy(RefinementFunctor* functor){refinement_functor_ = functor;};
	virtual ~RefinementStrategy(){};

	virtual void refine(GridStorage* storage, HashRefinementAbstract* hash_refinement)=0;

protected:

	RefinementFunctor* get_refinement_functor(){return refinement_functor_;}
	void set_refinement_functor(RefinementFunctor* functor) {refinement_functor_ = functor;}

private:
	RefinementFunctor* refinement_functor_;
};

}
}

#endif /* REFINEMENTSTRATEGY_HPP_ */
