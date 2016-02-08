/*
 * SerialCombiGrid.hpp
 *
 *  Created on: May 22, 2014
 *      Author: petzko
 */
#ifndef SERIALCOMBIGRID_HPP_
#define SERIALCOMBIGRID_HPP_

#include <sgpp/combigrid/combigrid/CombiGrid.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>
#include "../utils/combigrid_utils.hpp"

namespace combigrid {
template<typename _Tp>
class SerialCombiGrid: public CombiGrid<_Tp> {
 public:

  SerialCombiGrid(int in_dim,
                  const std::vector<bool>& in_hasBoundaryPts) : CombiGrid<_Tp>(in_dim,
                        in_hasBoundaryPts) {
    /**QUE NADA*/
  }

  /** sets the domain of all the full grids, this is the correct way for extrapolation */
  void setDomainAllFG(GridDomain* gridDomain) const;

  /** create the actual vector for the full grids. <br>
   * This is be different for serial and parallel implementations */
  void createFullGrids();


};

}

#endif /* SERIALCOMBIGRID_HPP_ */
