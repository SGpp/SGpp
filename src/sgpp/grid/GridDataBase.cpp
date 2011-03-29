/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "grid/GridDataBase.hpp"

#include "exception/file_exception.hpp"
#include "exception/data_exception.hpp"

#include "grid/GridStorage.hpp"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

namespace sg
{

  GridDataBase::GridDataBase(size_t dim) : _map(), _dim(dim)
  {}

  GridDataBase::GridDataBase(Grid *grid, DataVector &values) : _map(), _dim(grid->getStorage()->dim()) 
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
    grid_map_const_iterator git;

    for (git=_map.begin(); git != _map.end(); git++) {
      git->first->toString(ostream);
      ostream << " -> " << git->second << std::endl;
    }
    return ostream.str();
  }

  bool GridDataBase::hasKey(GridIndex* gi) {
    return _map.find(gi) != _map.end();
  }

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

  double GridDataBase::get(GridIndex* gi) {
    grid_map_const_iterator ind = _map.find(gi);
    if (ind == _map.end()) {
      //      throw new sg::data_exception("GridDataBase::get : grid point not in database");
      return NULL;
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

  void GridDataBase::save(std::string filename, bool type) {
    // streams
    std::ofstream fout;

    if (type==ascii) {

    // open file
    fout.open(filename.c_str());
    if (! fout.is_open()) {  
      std::string msg = "Error! Unable to open file '" + filename + "' for read access.";
      throw new file_exception(msg.c_str());
    }
    // dump contents to file
    // write ASCII information
    fout << type << std::endl;
    // write dimensionality
    fout << _dim << std::endl;
    // iterate over hashmap
    grid_map_const_iterator git;
    level_t lev; index_t ind;
    for (git=_map.begin(); git != _map.end(); git++) {
      for (int d=0; d<_dim; d++) {
	git->first->get(d, lev, ind);
	fout << lev << " " << ind << " ";
      }
      fout << git->second << std::endl;
    }

    } 
    else { // type: binary
    // open file
      fout.open(filename.c_str(), std::ios::binary);
    if (! fout.is_open()) {  
      std::string msg = "Error! Unable to open file '" + filename + "' for read access.";
      throw new file_exception(msg.c_str());
    }
    // dump contents to file
    // write binary information
    fout.write((char *)(&type), sizeof(type));
    // write dimensionality
    fout.write((char *)(&_dim), sizeof(_dim));
    // iterate over hashmap
    grid_map_const_iterator git;
    level_t lev; index_t ind;
    for (git=_map.begin(); git != _map.end(); git++) {
      for (int d=0; d<_dim; d++) {
	git->first->get(d, lev, ind);
	fout.write((char *)(&lev), sizeof(lev));
	fout.write((char *)(&ind), sizeof(ind));
      }
      fout.write((char *)&(git->second), sizeof(git->second));
    }


    }
    // tidy up
    fout.close();

  }

  void GridDataBase::load(const std::string filename) {
    // stream
    std::ifstream fin;
    bool ftype;
    int dim;
    _loadTypeDim(filename, ftype, dim, fin);

    std::cout << "ftype" << ftype << " dim " << dim << std::endl;
    
/*
    // open file
    fin.open(filename.c_str());
    if (! fin.is_open()) {  
      std::string msg = "GridDataBase::load Error! Unable to open file '" + filename + "' for read access.";
      throw new file_exception(msg.c_str());
    }
    // dimensionality
    int dim;
    fin >> dim;
    if (dim != _dim) {
      std::string msg = "GridDataBase::load Error! Dimensions do not match";
      throw new file_exception(msg.c_str());
    }
    // grid points
    GridIndex gi(_dim);
    level_t lev;
    index_t ind;
    double val;
    while (! fin.eof()) {
      for (int d=0; d<_dim; d++) {
	fin >> lev;
	fin >> ind;
	gi.set(d, lev, ind);
      }
      fin >> val;
      set(&gi, val);
    }
*/
  }

  void GridDataBase::_loadTypeDim(const std::string filename, bool &ftype, int &dim, std::ifstream &fin) {
    // open file, test in binary mode
    fin.open(filename.c_str(), std::ios::binary);
    if (! fin.is_open()) {  
      std::string msg = "GridDataBase: Error! Unable to open file '" + filename + "' for read access.";
      throw new file_exception(msg.c_str());
    }
    fin.read((char *)&ftype, sizeof(ftype));
    // we're in binary mode
    if (ftype == binary) {
      std::cout << filename << " is binary" << std::endl;
      fin.read((char *)&dim, sizeof(dim));
    }
    // we're in ASCII mode
    else {
      fin.close();
      fin.open(filename.c_str());
      std::cout << filename << " is ascii" << std::endl;
      fin >> ftype;
      if (ftype != ascii) {
	std::string msg = "GridDataBase: Error! Unrecognizable type of file '" + filename + "'";
	throw new file_exception(msg.c_str());
      }
      fin >> dim;
    }

/*

    // dimensionality
    int dim;
    fin >> dim;
    if (dim != _dim) {
      std::string msg = "GridDataBase::load Error! Dimensions do not match";
      throw new file_exception(msg.c_str());
    }
    // grid points
    GridIndex gi(_dim);
    level_t lev;
    index_t ind;
    double val;
    while (! fin.eof()) {
      for (int d=0; d<_dim; d++) {
	fin >> lev;
	fin >> ind;
	gi.set(d, lev, ind);
      }
      fin >> val;
      set(&gi, val);
    }
*/    



  }


}

