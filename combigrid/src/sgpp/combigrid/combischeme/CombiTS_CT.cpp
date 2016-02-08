#include <sgpp/combigrid/combischeme/CombiTS_CT.hpp>

template<typename _Tp>
combigrid::CombiTS_CT<_Tp>::CombiTS_CT(std::vector<int> minlevels,
		std::vector<int> maxlevels) {
	_levels_small = minlevels;
	this->_levels = maxlevels;
	this->_makeCombiInDimension.clear();
}

template<typename _Tp>
combigrid::CombiTS_CT<_Tp>::CombiTS_CT(const std::vector<int>& in_levels,
		const std::vector<bool>& makeCombiInDimension) {

	// set the makeCombiIndimensions flag vector
	this->_makeCombiInDimension = makeCombiInDimension;
	// set the small levels vector to be 0.5*in_levels_vector
	_levels_small.clear();
	for (unsigned int i = 0; i < in_levels.size(); i++)
		this->_levels_small.push_back(in_levels[i] / 2);

	// set the big (max) levels vector to be equal to in_levels
	this->_levels = in_levels;
}

template<typename _Tp>
combigrid::CombiTS_CT<_Tp>::CombiTS_CT(std::vector<int> in_max_levels) {
	// set the small levels vector to be 0.5*levels
	_levels_small.clear();
	for (unsigned int i = 0; i < in_max_levels.size(); i++)
		this->_levels_small.push_back(in_max_levels[i] / 2);
	// set the big (max) levels vector to be equal to in_levels
	this->_levels = in_max_levels;
	this->_makeCombiInDimension.clear();
}

template<typename _Tp>
void combigrid::CombiTS_CT<_Tp>::initCombiGrid(int in_dim,
		std::vector<std::vector<int> > & out_levels_vector,
		std::vector<_Tp>& out_coefs) {

	//initiate a duplicate
	//of the vectors with max and min levels values
	std::vector<int> levels_big = this->_levels;
	std::vector<int> levels_small = this->_levels_small;

	/* initiate the flag vector to indicate which are the active dimensions..
	 * counter of the number of active dimension! A dimension is considered active if (conditions are
	 * enumerated in order of decreasing strength if the ith condition is the first condition
	 * from the list to be satisfied, it invalidades all conditions below it !
	 *
	 * 1) the user has EXPLICITLY specified it as active (via the makeCombiInDimension array)
	 * 2) the user has specified maxlevels and minlevels array! then the current dimension is considered
	 * active if the maxlevel != minlevel for this dimension
	 * 3) default behavior --> all levels are active to begin with!
	 */

	std::vector<bool> active_dim(in_dim, true);
	int nr_activedimensions = in_dim;
	if (_makeCombiInDimension.size() != 0) {
		active_dim = _makeCombiInDimension;
	} else {
		if (levels_small.size() > 0 && levels_big.size() > 0)
			for (int i = 0; i < in_dim; i++)
				active_dim[i] = (levels_small[i] != levels_big[i]);
	}

	std::vector<int> level_tmp;

	// add all the spaces for all dimensions
	for (int i = 0; i < in_dim; i++) {
		if (active_dim[i]) {
			level_tmp = levels_small;
			level_tmp[i] = levels_big[i];
			out_levels_vector.push_back(level_tmp);
			out_coefs.push_back(1.0);
		} else {
			nr_activedimensions--;
		}
	}

	// add the smallest space
	out_levels_vector.push_back(levels_small);
	out_coefs.push_back(_Tp(1.0 - nr_activedimensions));
	this->removeDuplicates(out_levels_vector, out_coefs);
}

template<typename _Tp>
void combigrid::CombiTS_CT<_Tp>::re_initCombiGrid(int in_dim,
		const std::vector<FGridContainer<_Tp>*> in_grids,
		std::vector<std::vector<int> > & out_levels_vector,
		std::vector<_Tp>& out_coefs) {

	std::vector<std::vector<int> > tmp_levels_vector;
	std::vector<_Tp> tmp_coefs;
	initCombiGrid(in_dim, tmp_levels_vector, tmp_coefs);

	//examine what vectors do we have in the in_grids vector....
	//1) first deactivate all existing grids

	for (unsigned int i = 0; i < in_grids.size(); i++)
		in_grids[i]->deactivateGrid();

	//2) if the scheme attempts to create an already existing grid, change the coeffs and activate it.
	for (unsigned int i = 0; i < in_grids.size(); i++) {
		unsigned int j = 0;
		while (j < tmp_levels_vector.size()) {
			//if the grid vector i's levels vector == temp_levels_vector[j] and the coeffs are the same
			if (in_grids[i]->getFGLevels() == tmp_levels_vector[j])

			{	// if true leave the grid as activated..change its coefficient
				// and remove the record from the list of grids to be created.
				in_grids[i]->setCoef(tmp_coefs[j]);
				in_grids[i]->activateGrid();
				unsigned int nr = (unsigned int) tmp_levels_vector.size();
				tmp_levels_vector[j] = tmp_levels_vector[nr - 1];
				tmp_levels_vector.resize(nr - 1);
				tmp_coefs[j] = tmp_coefs[nr - 1];
				tmp_coefs.resize(nr - 1);
			} else {
				j++;
			}
		}
	}

	out_levels_vector = tmp_levels_vector;
	out_coefs = tmp_coefs;
	//now initialize the combigrid as if the in_grids vector were empty..

}

/** Implement this method with desired logic to handle situations when recomputation of the coefficients might be necessary.
 * 	Exemplary use cases could be the addition or removal of a fullgrid to/from the combigrid container.
 * @param in_dim - dimension of the problem
 * @param in_levels_vector - the vector containing the selected levels for each fullgrid of the caller combigrid
 * @param out_coeffs - (out) vector containing the combigrid coefficients (usually +1 or -1) -> the method stors the new
 * recomputed coefficients in the out_coeffs vector
 * */
template<typename _Tp>
void combigrid::CombiTS_CT<_Tp>::recomputeCoefficients(int in_dim,
		std::vector<FGridContainer<_Tp>*>& out_fgrids) {
	std::cout
			<< " combiTS_CT scheme -> recomputeCoefficients has been invoked \n";
} //action

template class combigrid::CombiTS_CT<float>;
template class combigrid::CombiTS_CT<double>;

// add more declarations at your hearth's will !!!!
