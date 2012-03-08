/*
 * RefinementStrategy.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: khakhutv_local
 */

#include "RefinementStrategy.hpp"
namespace sg
{
namespace base
{


RefinementStrategy::RefinementStrategy()
{
	// TODO Auto-generated constructor stub

}

RefinementStrategy::RefinementStrategy(RefinementFunctor* functor)
{
	refinement_functor_ = functor;
}



RefinementStrategy::~RefinementStrategy()
{
	// TODO Auto-generated destructor stub
}


}
}
