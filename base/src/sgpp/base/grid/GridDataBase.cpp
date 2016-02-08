// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridDataBase.hpp>

#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/base/grid/GridStorage.hpp>

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

GridDataBase::GridDataBase(size_t dim) : _map(), _dim(static_cast<int>(dim)) {}

GridDataBase::GridDataBase(Grid* grid, DataVector& values) : _map(),
  _dim(static_cast<int>(grid->getStorage()->dim())) {
  GridStorage* gs = grid->getStorage();

  for (size_t i = 0; i < gs->size(); i++) {
    set(gs->get(i), values[gs->seq(gs->get(i))]);
  }
}

GridDataBase::GridDataBase(const std::string& filename) {
  _map = grid_map();

  // stream
  std::ifstream fin;
  char ftype;
  int dim;
  // load
  _loadTypeDim(filename, ftype, dim, fin);
  _dim = dim;
  _loadData(fin, ftype);
}

GridDataBase::~GridDataBase() {
  clear();
}

void GridDataBase::clear() {
  for (grid_map_iterator ind = _map.begin(); ind != _map.end(); ind++) {
    // free allocated memory for GridIndex
    delete ind->first;
  }

  _map.clear();
}

std::string GridDataBase::toString() {
  // output stream
  std::ostringstream ostream;
  // iterator over hashmap
  grid_map_const_iterator git;

  for (git = _map.begin(); git != _map.end(); git++) {
    git->first->toString(ostream);
    ostream << " -> " << git->second << std::endl;
  }

  return ostream.str();
}


bool GridDataBase::hasKey(GridIndex* gi) {
  return _map.find(gi) != _map.end();
}

void GridDataBase::set(GridIndex* gi, float_t value) {
  grid_map_const_iterator ind = _map.find(gi);

  if (ind == _map.end()) {
    // create copy and store
    index_pointer cpy = new GridIndex(*gi);
    _map[cpy] = value;
  } else {
    // just overwrite value
    _map[ind->first] = value;
  }
}

void GridDataBase::setValuesFor(Grid* grid, DataVector& values) {
  GridStorage* gs = grid->getStorage();

  for (GridStorage::grid_map_iterator iter = gs->begin(); iter != gs->end();
       iter++) {
    values[iter->second] = get(iter->first);
  }
}

float_t GridDataBase::get(GridIndex* gi) {
  grid_map_const_iterator ind = _map.find(gi);

  if (ind == _map.end()) {
    std::cerr << gi->toString() << " not in database" << std::endl;
    throw new SGPP::base::data_exception("GridDataBase::get : grid point not in database");
  }

  return ind->second;
}

void GridDataBase::remove(GridIndex* gi) {
  grid_map_iterator ind = _map.find(gi);

  if (ind != _map.end()) {
    // remove
    _map.erase(ind);
    // free allocated memory for GridIndex
    delete ind->first;
  }
}

void GridDataBase::save(std::string filename, char ftype) {
  // streams
  std::ofstream fout;

  if (ftype == ascii) {
    // set precision
    fout << std::scientific;
    fout.precision(14);

    // open file
    fout.open(filename.c_str());

    if (!fout.is_open()) {
      std::string msg = "Error! Unable to open file '" + filename +
                        "' for read access.";
      throw new file_exception(msg.c_str());
    }

    // dump contents to file
    // write ASCII information
    fout << ftype << std::endl;
    // write dimensionality
    fout << _dim << std::endl;
    // iterate over hashmap
    grid_map_const_iterator git;
    level_t lev;
    index_t ind;

    for (git = _map.begin(); git != _map.end(); git++) {
      for (int d = 0; d < _dim; d++) {
        git->first->get(d, lev, ind);
        fout << lev << " " << ind << " ";
      }

      fout << git->second << std::endl;
    }

  } else { // type: binary
    // open file
    fout.open(filename.c_str(), std::ios::binary);

    if (!fout.is_open()) {
      std::string msg = "Error! Unable to open file '" + filename +
                        "' for read access.";
      throw new file_exception(msg.c_str());
    }

    // dump contents to file
    // write binary information
    fout.write((char*)(&ftype), sizeof(ftype));
    // write dimensionality
    fout.write((char*)(&_dim), sizeof(_dim));
    // iterate over hashmap
    grid_map_const_iterator git;
    level_t lev;
    index_t ind;

    for (git = _map.begin(); git != _map.end(); git++) {
      for (int d = 0; d < _dim; d++) {
        git->first->get(d, lev, ind);
        fout.write((char*)(&lev), sizeof(lev));
        fout.write((char*)(&ind), sizeof(ind));
      }

      fout.write((const char*) & (git->second), sizeof(git->second));
    }
  }

  // tidy up
  fout.close();
}

void GridDataBase::load(const std::string filename) {
  // stream
  std::ifstream fin;
  char ftype;
  int dim;
  _loadTypeDim(filename, ftype, dim, fin);
  std::cout << "ftype" << ftype << " dim " << dim << std::endl;

  if (dim != _dim) {
    std::string msg = "GridDataBase::load Error! Dimensions do not match";
    throw new file_exception(msg.c_str());
  }

  _loadData(fin, ftype);
}

void GridDataBase::_loadTypeDim(const std::string filename, char& ftype,
                                int& dim, std::ifstream& fin) {
  // open file, test in binary mode
  fin.open(filename.c_str(), std::ios::binary);

  if (!fin.is_open()) {
    std::string msg = "GridDataBase: Error! Unable to open file '" + filename +
                      "' for read access.";
    throw new file_exception(msg.c_str());
  }

  fin.read((char*)&ftype, sizeof(ftype));

  // we're in binary mode
  if (ftype == binary) {
    std::cout << filename << " is binary" << std::endl;
    fin.read((char*)&dim, sizeof(dim));
  }
  // we're in ASCII mode
  else {
    fin.close();
    fin.open(filename.c_str());
    std::cout << filename << " is ascii" << std::endl;
    fin >> ftype;

    if (ftype != ascii) {
      std::string msg = "GridDataBase: Error! Unrecognizable type of file '" +
                        filename + "'";
      throw new file_exception(msg.c_str());
    }

    fin >> dim;
  }

}

void GridDataBase::_loadData(std::ifstream& fin, char& ftype) {
  // check if file handler is allright
  if (!fin.is_open()) {
    return;
  }

  // read grid points
  GridIndex gi(_dim);
  level_t lev;
  index_t ind;
  float_t val;

  // ASCII data
  if (ftype == ascii) {
    while (!fin.eof()) {
      for (int d = 0; d < _dim; d++) {
        fin >> lev;
        fin >> ind;
        gi.set(d, lev, ind);
      }

      fin >> val;
      set(&gi, val);
    }
  }
  // binary data
  else {
    while (!fin.eof()) {
      for (int d = 0; d < _dim; d++) {
        fin.read((char*)&lev, sizeof(lev));
        fin.read((char*)&ind, sizeof(ind));
        gi.set(d, lev, ind);
      }

      fin.read((char*)&val, sizeof(val));
      set(&gi, val);
    }
  }

  // close fin
  fin.close();
}

GridDataBase::grid_map_iterator GridDataBase::begin() {
  return _map.begin();
}

GridDataBase::grid_map_iterator GridDataBase::end() {
  return _map.end();
}

}  // namespace base
}  // namespace SGPP
