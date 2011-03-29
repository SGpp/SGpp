/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef GRIDDATABASE_HPP
#define GRIDDATABASE_HPP

#include "common/hash_map_config.hpp"

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"

#include <string>

namespace sg
{

/**
 * Class that allows to keep a storage of arbitrary grid points. It has
 * the functionality of a dictionary that is used for storing and 
 * retrieving grid points. Internally, a hash_map is used.
 *
 * Note: GridDataBase currently supports only to store pairs of
 * (grid point -> double).
 */
class GridDataBase
{
private:
  typedef GridIndex index_type;
  typedef GridIndex* index_pointer;
#ifndef USETRONE
#ifndef LARRABEENATIVE
  typedef std::hash_map<index_pointer, double, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
#endif
#ifdef LARRABEENATIVE
  typedef std::hash_map<index_pointer, double, LRBSGHasher<index_pointer> > grid_map;
#endif
#endif
#ifdef USETRONE
  typedef std::tr1::unordered_map<index_pointer, double, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
#endif

  // the hash_map
  grid_map _map;
  // dimensionality of grid
  int _dim;

  // index and level types
  typedef GridIndex::index_type index_t;
  typedef GridIndex::level_type level_t;

  /**
   * Loads database in ASCII format from file. Adds (grid point - value) mappings
   * to current database. Overwrites existing entries. To load a new database, 
   * use GridDataBase::fromFile(std::string filename).
   * @param filename name of file
   */
  void _loadTypeDim(const std::string filename, bool &ftype, int &dim, std::ifstream &fin);

  /**
   * Loads database in ASCII format from file. Adds (grid point - value) mappings
   * to current database. Overwrites existing entries. To load a new database, 
   * use GridDataBase::fromFile(std::string filename).
   * @param filename name of file
   */
  void _loadGrid(std::ifstream fin);



public:

  typedef grid_map::iterator grid_map_iterator;
  typedef grid_map::const_iterator grid_map_const_iterator;

  static const bool ascii=true;
  static const bool binary=false;

  /**
   * Standard Constructor, creating an empty database with dimensionality dim.
   * @param dim the dimensionality of the grid points
   */
  GridDataBase(size_t dim);

  /**
   * Constructor, copying from existing grid and values.
   * @param grid the grid to copy from
   * @param values the initial values
   */
  GridDataBase(Grid *grid, DataVector &values);

  /**
   * Constructor, reading grid data base from file.
   * @param filename name of file
   */
  //  GridDataBase(std::string filename);

  /**
   * Returns std::string representation of database.
   * @return string representation
   */
  std::string toString();

  /**
   * Test, whether database already contains a grid point.
   * @param gi a grid index
   * @return false, if grid point not in database
   */
  bool hasKey(GridIndex* gi);

  /**
   * Store grid point - value pair in database.
   * @param gi a grid index
   * @param value the value to store
   */
  void set(GridIndex* gi, double value);

    /**
     * Returns the number of grid points that are stored
     * in the database.
     *
     * @return the size of the database
     */
    size_t size() const
    {
        return _map.size();
    };

    /**
     * Returns the dimensionality of the grid points.
     *
     * @return the dimensionality of the grid points
     */
    size_t dim() const
    {
      return _dim;
    };

  /**
   * Get value for grid point from database.
   * @param gi a grid index
   * @return value
   */
  double get(GridIndex* gi);

  /**
   * Remove grid point from database. Do nothing, if not in database.
   * @param gi grid point
   */
  void remove(GridIndex* gi);

  /**
   * Save database in ASCII format to file. Overwrites existing files
   * without warning! Each line contains an entry @f$l_1, i_1, \ldots, l_d, i_d, value@f$.
   * @param filename name of file
   * @param type filetype (either ASCII, GridDataBase::asc (default), or binary, GridDataBase::bin)
   */
  void save(std::string filename, bool type=GridDataBase::ascii);

  /**
   * Loads database in ASCII format from file. Adds (grid point - value) mappings
   * to current database. Overwrites existing entries. To load a new database, 
   * use GridDataBase::fromFile(std::string filename).
   * @param filename name of file
   */
  void load(const std::string filename);


};

}
#endif /* GRIDDATABASE_HPP */

