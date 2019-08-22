// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRIDDATABASE_HPP
#define GRIDDATABASE_HPP

// #include <sgpp/base/tools/hash_map_config.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <string>
#include <unordered_map>


namespace sgpp {
namespace base {

/**
 * Class that allows to keep a storage of arbitrary grid points. It has
 * the functionality of a dictionary that is used for storing and
 * retrieving grid points. Internally, a hash_map is used.
 *
 * Note: GridDataBase currently supports only to store pairs of
 * (grid point -> double).
 */
class GridDataBase {
 private:
  typedef GridPoint index_type;
  typedef GridPoint* index_pointer;
  /*
  #ifndef USETRONE
  #ifndef LARRABEENATIVE
  typedef std::hash_map<index_pointer, double, sgpp::base::hash<index_pointer>, sgpp::base::eqIndex<index_pointer> > grid_map;
  #endif
  #ifdef LARRABEENATIVE
  typedef std::hash_map<index_pointer, double, LRBSGHasher<index_pointer> > grid_map;
  #endif
  #endif
  #ifdef USETRONE
  typedef std::tr1::unordered_map<index_pointer, double, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
  #endif
  */
  //        typedef HashGridStorage::grid_map grid_map;
  typedef std::unordered_map<index_pointer, double,
          HashGridPointPointerHashFunctor,
          HashGridPointPointerEqualityFunctor >
          grid_map;

  // the hash_map
  grid_map _map;
  // dimensionality of grid
  int _dim;

  // index and level types
  typedef GridPoint::index_type index_t;
  typedef GridPoint::level_type level_t;

  /**
   * Loads database in ASCII format from file. Adds (grid point - value) mappings
   * to current database. Overwrites existing entries. To load a new database,
   * use GridDataBase::fromFile(std::string filename).
   * @param filename name of file
   * @param ftype type of stream
   * @param dim dimension
   * @param fin file stream
   */
  void _loadTypeDim(const std::string filename, char& ftype, int& dim,
                    std::ifstream& fin);

  /**
   * Loads database in ASCII format from file. Adds (grid point - value) mappings
   * to current database. Overwrites existing entries. To load a new database,
   * use GridDataBase::fromFile(std::string filename).
   * closes filestream fin at end
   * @param fin file stream
   * @param ftype type of stream
   */
  void _loadData(std::ifstream& fin, char& ftype);

 public:
  typedef grid_map::iterator grid_map_iterator;
  typedef grid_map::const_iterator grid_map_const_iterator;

  static const char ascii = 'a';
  static const char binary = 'b';

  /**
   * Standard Constructor, creating an empty database with dimensionality dim.
   * @param dim the dimensionality of the grid points
   */
  explicit GridDataBase(size_t dim);

  /**
   * Constructor, copying from existing grid and values.
   * @param grid the grid to copy from
   * @param values the initial values
   */
  GridDataBase(Grid* grid, DataVector& values);

  /**
   * Constructor, reading from existing database.
   * @param filename filename of database file
   */
  explicit GridDataBase(const std::string& filename);

  /**
   * Destructor
   */
  ~GridDataBase();

  /**
   * Clears database, removing all entries. Dimensionality is maintained.
   */
  void clear();

  /**
   * Returns std::string representation of database.
   * @return string representation
   */
  std::string toString();

  /**
   * Test, whether database already contains a grid point.
   * @param gi a grid point
   * @return false, if grid point not in database
   */
  bool hasKey(GridPoint* gi);

  /**
   * Store (grid point - value) pair in database.
   * @param gi a grid point
   * @param value the value to store
   */
  void set(GridPoint* gi, double value);

  /**
   * Set values for given grid from stored ones in database.
   * @param grid a grid
   * @param values the corresponding coefficient vector
   */
  void setValuesFor(Grid* grid, DataVector& values);

  /**
   * Returns the number of grid points that are currently stored
   * in the database.
   *
   * @return the size of the database
   */
  size_t size() const {
    return _map.size();
  }

  /**
   * Returns the dimensionality of the grid points.
   *
   * @return the dimensionality of the grid points
   */
  size_t dim() const {
    return _dim;
  }

  /**
   * Get value for grid point from database. Return nullptr if not existant.
   * @param gi a grid point
   * @return value
   */
  double get(GridPoint* gi);

  /**
   * Remove grid point from database. Do nothing, if not in database.
   * @param gi grid point
   */
  void remove(GridPoint* gi);

  /**
   * Save database to file. Overwrites existing files
   * without warning! Writes either ASCII (default) or binary format.
   * @param filename name of file
   * @param ftype filetype (either ASCII, GridDataBase::asc (default), or binary, GridDataBase::bin)
   */
  void save(std::string filename, char ftype = GridDataBase::ascii);

  /**
   * Loads database in ASCII or binary format. Adds (grid point - value) mappings
   * to current database. Overwrites existing entries. To load a new database,
   * use GridDataBase::GridDataBase(std::string filename).
   * @param filename name of file
   */
  void load(const std::string filename);

  /**
   * Return iterator to beginning. Entries are of type pair<GridPoint, double>.
   * @return iterator to beginning.
   */
  grid_map_iterator begin();

  /**
   * Return iterator to end.
   * @return iterator to end.
   */
  grid_map_iterator end();
};

}  // namespace base
}  // namespace sgpp
#endif /* GRIDDATABASE_HPP */
