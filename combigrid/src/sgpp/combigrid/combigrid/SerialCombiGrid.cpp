

#include <sgpp/combigrid/combigrid/SerialCombiGrid.hpp>

/** sets the domain of all the full grids, this is the correct way for extrapolation */
template <typename _Tp>
void combigrid::SerialCombiGrid< _Tp >::setDomainAllFG(GridDomain* gridDomain) const{
	for (unsigned int i = 0; i < this->_fgrids.size(); i++) {
		this->_fgrids[i]->fg()->setDomain(gridDomain);
	}
}

/** create the actual vector for the full grids. <br>
 * This is be different for serial and parallel implementations */
template <typename _Tp>
void combigrid::SerialCombiGrid<_Tp>::createFullGrids() {

	for (unsigned int i = 0; i < this->_fgrids.size(); i++) {
		this->_fgrids[i]->createFullGrid();
	}
}


template class combigrid::SerialCombiGrid<float>;
template class combigrid::SerialCombiGrid<double>;

// add more declarations at your hearth's will !!!!
