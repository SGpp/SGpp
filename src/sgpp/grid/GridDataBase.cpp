/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "grid/GridDataBase.hpp"

#include "exception/data_exception.hpp"

#include "grid/GridStorage.hpp"

#include <string>
#include <sstream>

namespace sg
{

  GridDataBase::GridDataBase(size_t dim) : _dim(dim), _map() 
  {}

  GridDataBase::GridDataBase(Grid *grid, DataVector &values) : _dim(grid->getStorage()->dim()), _map() 
  {
    GridStorage* gs = grid->getStorage();
    for (size_t i = 0; i < gs->size(); i++) {
      set(gs->get(i), values[gs->seq(gs->get(i))]);
    }
  }

  std::string GridDataBase::toString() {
    // output stream
    std::ostringstream ostream;
    // iterator over hashmap
    grid_map_iterator git;

    for (git=_map.begin(); git != _map.end(); git++) {
      git->first->toString(ostream);
      ostream << " -> " << git->second << std::endl;
    }
    return ostream.str();
  }

  /**
   * Test, whether database already contains a grid point.
   * @param gi a grid index
   * @return false, if grid point not in database
   */
  bool GridDataBase::hasKey(GridIndex* gi) {
    return _map.find(gi) != _map.end();
  }

  /**
   * Store grid point - value pair in database.
   * @param gi a grid index
   * @param value the value to store
   */
  void GridDataBase::set(GridIndex* gi, double value) {
    grid_map_const_iterator ind = _map.find(gi);
    if (ind == _map.end()) {
      // create copy and store
      index_pointer cpy = new GridIndex(gi);
      _map[cpy] = value;
    } else {
      // just overwrite value
      _map[ind->first] = value;
    }
  }


  /**
   * Get value for grid point from database.
   * @param gi a grid index
   * @return value
   */
  double GridDataBase::get(GridIndex* gi) {
    grid_map_const_iterator ind = _map.find(gi);
    if (ind == _map.end()) 
      throw new sg::data_exception("GridDataBase::get : grid point not in database");
    return ind->second;
  }

  void GridDataBase::remove(GridIndex* gi) {
    //    grid_map_const_iterator ind = _map.find(gi);
    //    if (ind != _map.end()) 
    _map.erase(gi);
  }



}

