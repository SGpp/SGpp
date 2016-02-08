/*
 * CombiArbitraryScheme_demo.hpp
 *
 *  Created on: 29 Jun 2014
 *      Author: kenny
 */

#ifndef COMBIARBITRARYSCHEME_HPP_
#define COMBIARBITRARYSCHEME_HPP_

#include <sgpp/combigrid/combischeme/AbstractCombiScheme.hpp>
#include <sgpp/combigrid/utils/CombigridLevelVector.hpp>

namespace combigrid {
/** Combischeme created with from an arbitrary active set of full grids. Missing
 * subgrids for the creation of a valid combi solution are automatically added.
 */
template<typename _Tp>
class CombiArbitraryScheme: public combigrid::AbstractCombiScheme<_Tp> {

private:
	std::vector<std::vector<int> > _levels_vector;

public:

	CombiArbitraryScheme(std::vector< std::vector<int> > levels_vector);

	/**
	 * do nothing constructor
	 */
	virtual ~CombiArbitraryScheme() {
		;
	}

	/** The interface function that all classes inheriting from AbstractCombiScheme should implement.
	 * initCombiGrid implements the desired combiGrid scheme, i.e. method of construction of the functional spaces
	 * and the corresponding coefficient for each space. Notice that the combischeme DOES NOT own any of the data
	 * that is given to it. It simply fills in the "out_levels_vector" and "out_coefs" vector data and returns without
	 * keeping an internal copies of it. This is a nice way to enforce decoupling between application logic implementation and
	 * data management!!
	 *
	 *
	 * @param in_dim (in) an input parameter specifying the dimension of the problem
	 * @param out_levels_vector (out) a vector containing the selected hierarchical levels by the current combigrid implementation
	 * @param out_coefs (out) a vector containing the interpolation coeficients for each level...
	 *
	 * */

	void initCombiGrid(int in_dim,
			std::vector<std::vector<int> > & out_levels_vector,
			std::vector<_Tp>& out_coefs);

	/**
	 *	re_initCombigrid - a method similar to the initCombiGrid, but this one allows data reusability.
	 *	The combischeme will examine the already existing full grids and as a result decide which of the already
	 *	existing grids to activate, which to deactivate and which to create anew! The output out_levels_vector and out_coefs vectors
	 *	contain the information ONLY for the ADDITIONAL grids that need to be created.
	 *
	 * @param in_dim (in) an input parameter specifying the dimension of the problem
	 * @param in_grids the already existing fullgrids. The scheme uses the information contained inside "in_grids" to decide which grids to activate and which to deactivate
	 * @param out_levels_vector (out) a vector containing the selected levels by the current combigrid implementation
	 * @param out_coefs (out) a vector containing the interpolation coefficients for each level...
	 *
	 * */

	void re_initCombiGrid(int in_dim,
			const std::vector<FGridContainer<_Tp>*> in_grids,
			std::vector<std::vector<int> > & out_levels_vector,
			std::vector<_Tp>& out_coefs);

	/** Implement this method with desired logic to handle situations when recomputation of the coefficients might be necessary.
	 * 	Exemplary use cases could be the addition or removal of a fullgrid to/from the combigrid container.
	 * @param in_dim - dimension of the problem
	 * @param out_fgrids - a vector containing the fullgrids associated with the current combigrid- will be updated with the newly recomputed
	 * coefficients
	 * */
	void recomputeCoefficients(int in_dim,
			std::vector<FGridContainer<_Tp>*>& out_fgrids);


};

}

#endif /* COMBIARBITRARYSCHEME_DEMO_HPP_ */
