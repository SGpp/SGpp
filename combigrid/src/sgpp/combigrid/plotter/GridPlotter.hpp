// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRIDPLOTTER_HPP_
#define GRIDPLOTTER_HPP_

#include <sgpp/combigrid/combigrid/AbstractCombiGrid.hpp>

namespace combigrid {

class Evaluable {

 public:
  Evaluable(const FullGridD* fg): fg_(fg) , cg_(0) {
    ;
  }
  Evaluable(const AbstractCombiGrid* cg): fg_(0) , cg_(cg) {
    ;
  }

  inline double eval(std::vector<double>& coords ) const {
    if (fg_ == 0) {
      return cg_->eval(coords);
    } else {
      return fg_->eval(coords);
    }
  }
 private:
  const FullGridD* fg_;
  const AbstractCombiGrid* cg_;
};

/** plott one grid in a matlab file */
class GridPlotter {
 public:

  /** empty Ctor*/
  GridPlotter() {
    ;
  }

  /** empty Dtor*/
  virtual ~GridPlotter() {
    ;
  }

  /** plot one Full grid */
  static void plotFullGrid(const std::string& filePath , const FullGridD* fg ,
                           std::vector<double>& globalCoord_in , int resolution = 0);

  /** plot one combination grid */
  static void plotCombiGrid(const std::string& filePath ,
                            const AbstractCombiGrid* cg ,
                            std::vector<double>& globalCoord_in , int resolution = 0);

 private:

  static void plotObject(int dim ,
                         const std::string& filePath ,
                         const Evaluable* obj ,
                         const GridDomain* domain ,
                         std::vector<double>& globalCoord_in ,
                         int resolution);

};

}

#endif /* GRIDPLOTTER_HPP_ */