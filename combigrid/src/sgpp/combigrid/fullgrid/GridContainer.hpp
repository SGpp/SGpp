/*
 * GridContainer.hpp
 *
 *  Created on: 23 Jun 2014
 *      Author: kenny
 */

#ifndef GRIDCONTAINER_HPP_
#define GRIDCONTAINER_HPP_


#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>
#include "../utils/combigrid_utils.hpp"


/***
 *
 *
 * Simple class that encapsulates the combigrid relevant full grid data. Those are:
 * - the fullgrid datatype itself
 * - the combination technique coefficient
 * - the levels vector for the current fullgrid!
 *
 */
namespace combigrid {

template<typename _Tp>
class FGridContainer {

 private:

  /***
   * levels vector (specifying level depth for each dimension) for the corresponding grid!
   */
  std::vector<int> _fullgrid_levels;
  /***
   * Combigrid coefficient! Special _Tp type!
   *
   */
  _Tp _coefficient;
  /**
   * fullgrid container!
   */
  FullGrid<_Tp>* _fullgrid;

  /**
   *
   * flag determining whether or not this grid values hsould be included in the evaluation
   *
   */
  bool _active;

 public:

  void setCoef(_Tp newCoef) {
    _coefficient = newCoef;
  }

  /**
   * Returns the coefficient of the underlying fullgrid
   *
   */
  _Tp getCoef() {
    return _coefficient;
  }

  void setFGLevels(std::vector<int> levels) {
    _fullgrid_levels = levels;
  }


  /***
   * Returns the levels vector of this fullgrid
   */
  std::vector<int> getFGLevels() {
    return _fullgrid_levels;
  }

  int getMaxLevel() {


    int max = 1;
    int DIM = (int) _fullgrid_levels.size();

    for (int d = 0 ; d < DIM; d++)
      max =  _fullgrid_levels[d] > max ? _fullgrid_levels[d] : max;

    return max;
  }



  /**
   * Returns the underlying fullgrid
   *
   */
  FullGrid<_Tp>* fg() {
    return _fullgrid;
  }

  void activateGrid() {
    _active = true;
  }
  void deactivateGrid() {
    _active = false;
  }
  bool isActive() {
    return _active;
  }

  /**
   * FGridContainer constructor. No empty, default, constructors allowed!
   *
   * @param levels an integer list that specifies the hierarchical level depth for each dimension.
   * @param hasBoundaryPts : a boolean list with flags specifying if a particular dimension has boundary points included!
   * @param coef: well the combigrid coefficient
   * @param basis: the basis function for the underlying fullgrid! This parameter has a default value NULL, in which case
   * the default LinearBasisFunction is set!
   *
   */
  FGridContainer(const std::vector<int> levels,
                 const std::vector<bool>& hasBoundaryPts, _Tp coef,
                 BasisFunctionBasis* basis = NULL) {

    int dim = static_cast<int>(levels.size());

    _fullgrid = new combigrid::FullGrid<_Tp>(dim, levels, hasBoundaryPts, basis);
    _fullgrid_levels = levels;
    _coefficient = coef;
    _active = true;

  }

  /***
   * Free memory allocated for the fullgrid
   *
   *
   */
  ~FGridContainer() {
    _fullgrid_levels.clear();
    delete _fullgrid;
  }


  /**
   *
   * Allocate actual memory for the fullgrid elements! A simple call to Fullgrid->createFullGrid();
   *
   */
  void createFullGrid() {
    if (_fullgrid != NULL)
      _fullgrid->createFullGrid();
  }

};

}

#endif /* GRIDCONTAINER_HPP_ */

