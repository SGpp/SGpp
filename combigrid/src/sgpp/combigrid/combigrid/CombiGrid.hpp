/*
 * CombiGrid.hpp
 *
 *  Created on: May 22, 2014
 *      Author: petzko
 */

#ifndef COMBIGRID_HPP_
#define COMBIGRID_HPP_

#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>

// ------ SGpp includes -------------
//#include "base/datatypes/DataVector.hpp"
//#include "base/grid/GridStorage.hpp"
#include <sgpp/combigrid/combischeme/AbstractCombiScheme.hpp>
#include <sgpp/combigrid/fullgrid/GridContainer.hpp>
#include <stdlib.h>
#include "../utils/combigrid_utils.hpp"
#include "../utils/combigrid_utils.hpp"


namespace combigrid {

template<typename _Tp>
class CombiGrid {

 public:


  CombiGrid(int in_dim, const std::vector<bool>& in_hasBoundaryPts) {

    _dim = in_dim;
    _fgrids.resize(0);
    _hasBoundaryPts = in_hasBoundaryPts;
    _combischeme = NULL;
    _domainMax.resize(_dim, 1.0);
    _domainMin.resize(_dim, 0.0);
    _stretchingType = EQUIDISTANT;
  }

  /**
   * Dtor which deletes all data of
   * associated with this combigrid
   */
  virtual ~CombiGrid() {
    deleteFullGrids();
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //////////////////////// Public Data management methods ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////

  /** evaluate the combi grid at one specified point
   * @param coords the coordinates to evaluate
   * @return the result from the combigrid evaluation
   * */
  virtual _Tp eval(const std::vector<double>& coords) const {

    _Tp res = 0.0;

    if (_combischeme == NULL) {
      COMBIGRID_OUT_WRN(
        "Combigrid Evaluation failed! Combischeme not set.",
        __FILE__, __LINE__)

    } else {
      res = _combischeme->evalCombiGrid( _fgrids, coords);
    }

    return res;
  }

  /** evaluate the combi grid at one specified point. Buffered evaluation
   * which might be faster than one evaluation point.
   * @param coords - a list of coordinates to evaluate the grid on
   * @param results the result vector
   */
  virtual void eval(std::vector<std::vector<double> >& coords,
                    std::vector<_Tp>& results) const {

    if (_combischeme == NULL) {
      COMBIGRID_OUT_WRN(
        "Combigrid Evaluation failed! Combischeme not set.",
        __FILE__, __LINE__)
    } else {
      _combischeme->evalCombiGrid( _fgrids, coords,
                                   results);
    }

  }

  /**
   * delete all the fullgrids contained in this combigrid
   */
  void deleteFullGrids() {

    for (unsigned int i = 0; i < _fgrids.size(); i++) {
      delete _fgrids[i];
    }

    _fgrids.resize(0);
  }

  /** delete the i-th full grid
   * @param i
   *
   */
  void deleteFullGrid(unsigned int i) {

    //delete one full grid
    size_t nrFG = _fgrids.size();

    if (i < nrFG) {
      // here we swap the last element with the i-th element
      FGridContainer<_Tp>* fg = _fgrids[i];
      _fgrids[i] = _fgrids[nrFG - 1];
      _fgrids[nrFG - 1] = fg;
      delete _fgrids[nrFG - 1];
      _fgrids.resize(nrFG - 1);

      // on state change --> call the combischeme recompute coefficients method!
      if (_combischeme != NULL)
        COMBIGRID_OUT_WRN(
          "Combigrid coefficients recalculation failed! Combischeme not set.",
          __FILE__, __LINE__)
        else
          _combischeme->recomputeCoefficients(_dim, _fgrids);

    }
  }

  /** adds a full grid with the specified level (dimension was specified in the Ctor) and boundary
   * @param levels levels of the FG
   * @param hasBoundaryPts for each dimension if the full grid should have boundary points
   * @param coef the coeficient in the combination scheme */
  void addFullGrid(const std::vector<int>& levels,
                   const std::vector<bool>& hasBoundaryPts, _Tp coef) {

    combigrid::FGridContainer<_Tp>* grid = new FGridContainer<_Tp>(levels,
        hasBoundaryPts, coef);
    _fgrids.push_back(grid);

    // on state change --> call the combischeme recalculate coefficients method!
    if (_combischeme == NULL)
      COMBIGRID_OUT_WRN(
        "Combigrid coefficients recalculation failed! Combischeme not set.",
        __FILE__, __LINE__)
      else
        _combischeme->recomputeCoefficients(_dim, _fgrids);

  }

  /**
   *
   * method which initializes the
   * combigrid with a set of full grids and the corresponding coefficients
   * !! The actual computation of the hierarchical spaces and the combigrid coefficients
   * is delegated to the underlying combischeme object!!!
   *
   *
   */
  void init() {

    if (_combischeme == NULL) {
      COMBIGRID_OUT_WRN(
        "Combigrid Evaluation failed! Combischeme not set.",
        __FILE__, __LINE__)
      return;
    }

    // first delete everything what was there before
    deleteFullGrids();
    std::vector<std::vector<int> > tmp_levels_vtr;
    std::vector<_Tp> tmp_cf_vector;

    // then reinitialize the levels vector and the coefs vector;
    _combischeme->initCombiGrid(_dim, tmp_levels_vtr, tmp_cf_vector);

    // create add all full grids to the fullgrid container.
    for (unsigned int i = 0; i < tmp_levels_vtr.size(); i++) {
      FGridContainer<_Tp>* fg = new FGridContainer<_Tp>(tmp_levels_vtr[i],
          _hasBoundaryPts, tmp_cf_vector[i]);
      _fgrids.push_back(fg);

    }

  }

  /**
   *
   * !! The actual computation of the hierarchical spaces and the combigrid coefficients
   * is delegated to the underlying combischeme object!!!
   *
   *  re_init - a method similar to init(), but this one allows data reusability!!
   *  Again the computation of the spaces
   *  and coefficients is delegated to the underlying combischeme, which  will
   *  examine the already existing full grids and as a result decide which of them to
   *  activate, which to deactivate and which to create anew!
   *
   */

  void re_init() {

    if (_combischeme == NULL) {
      COMBIGRID_OUT_WRN(
        "Combigrid Evaluation failed! Combischeme not set.",
        __FILE__, __LINE__)
      return;
    }

    // first delete everything what was there before
    std::vector<std::vector<int> > tmp_levels_vtr;
    std::vector<_Tp> tmp_cf_vector;

    // then reinitialize the levels vector and the coefs vector;
    _combischeme->re_initCombiGrid(_dim, _fgrids, tmp_levels_vtr,
                                   tmp_cf_vector);

    // create add all full grids to the fullgrid container.
    for (unsigned int i = 0; i < tmp_levels_vtr.size(); i++) {
      FGridContainer<_Tp>* fg = new FGridContainer<_Tp>(tmp_levels_vtr[i],
          _hasBoundaryPts, tmp_cf_vector[i]);
      _fgrids.push_back(fg);

    }

  }

  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////Getters & Setters//////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  std::vector<_Tp> getCoefs() {
    std::vector<_Tp> vect;

    for (unsigned int i = 0; i < _fgrids.size(); i++)
      vect.push_back(_fgrids[i]->getCoef());

    return vect;
  }

  std::vector<std::vector<int> > getLevelsVector() {
    std::vector<std::vector<int> > vect;

    for (unsigned int i = 0; i < _fgrids.size(); i++)
      vect.push_back(_fgrids[i]->getFGLevels());

    return vect;
  }

  const std::vector<bool>& getBoundaryFlags() const {
    return _hasBoundaryPts;
  }

  /**
   * Iterates over all active grids of the combigrid and sets the grid domain as
   * a homogeneously stretched domain spawning the [min;max] dim-dimensional box
   *
   * @param min - a vector of length dim, which contains the left boundaries of the
   * domain in each direction (dimension)
     * @param max - a vector of length dim, which contains the right boundaries of the
   * domain in each direction (dimension)
   * @param stretching - a pointer to a stretching class instance
   *
   * */
  void initializeActiveGridsDomain(std::vector<double> min,
                                   std::vector<double> max,
                                   combigrid::AbstractStretchingMaker* stretching) {

    for (unsigned int g = 0; g < this->getNrFullGrids(); g++) {
      if (getFullGrid(g)->isActive()) {
        FullGrid<_Tp>* fg = getFullGrid(g)->fg();
        GridDomain* domain = fg->getDomain();

        if (domain != NULL)
          delete domain;

        fg->setDomain(
          new GridDomain(_dim, fg->getLevels(), min, max,
                         *stretching));
      }
    }

    _domainMax = max ;
    _domainMin = min;

    if (stretching != NULL)
      _stretchingType = stretching->getStretchingType();
    else
      _stretchingType = EQUIDISTANT;
  }
  /**
   * Iterates over all active grids of the combigrid and sets the grid domain to dom
   *
   * @param dom - the new GrdiDomain
   *
   */
  void initializeActiveGridsDomain(GridDomain* dom) {

    for (int g = 0; g < this->getNrFullGrids(); g++) {
      if (getFullGrid(g)->isActive()) {
        FullGrid<_Tp>* fg = getFullGrid(g)->fg();
        GridDomain* domain = fg->getDomain();

        if (domain != NULL)
          delete domain;

        fg->setDomain(dom);
      }
    }

    _domainMax.clear();
    _domainMin.clear();

    for (int d = 0 ; d < _dim ; d++) {
      Domain1D _1ddom = dom->get1DDomain(d);
      _domainMax.push_back(_1ddom.getMaxDomain());
      _domainMin.push_back(_1ddom.getMinDomain());

    }

    _stretchingType = dom->getStretchingType();
  }

  Stretching getStretchingType() {
    return _stretchingType;
  }

  /** return the combi scheme */
  const AbstractCombiScheme<_Tp>* getCombiScheme() const {
    return _combischeme;
  }

  /**
   * Attach a new combischeme to this grid! An combischeme captures the application logic
   * (evaluation, fullgrid initialization, coefficient computation etc.) while the combigrid serves
   * only as a data storage entity. All of the aforementioned operations that have to be performed on a combigrid,
   * are delegated to the underlying combischeme!
   *
   */
  void attachCombiScheme(AbstractCombiScheme<_Tp>* scheme) {
    this->_combischeme = scheme;
  }

  /*
   * Detach the current combischeme, and return it to the user!
   *
   */
  AbstractCombiScheme<_Tp>* detachCombiScheme(void) {
    AbstractCombiScheme<_Tp>* res = _combischeme;
    this->_combischeme = NULL;
    return res;

  }
  /** return the number of full grids */
  unsigned int getNrFullGrids() const {
    return (unsigned int)_fgrids.size();
  }
  /**
   * Count the number of grids set to active!
   */
  unsigned int getNrActiveFullGrids() {
    unsigned int ctr = 0;

    for (unsigned int i = 0; i < _fgrids.size(); i++)
      ctr += _fgrids[i]->isActive() ? 1 : 0;

    return ctr;
  }

  /** return the full grids level vector */
  std::vector<int> getFullGridLevel(unsigned int i) {

    if (i >= _fgrids.size())
      return std::vector<int>();
    else
      return _fgrids[i]->getFGLevels();

  }

  FGridContainer<_Tp>* getFullGrid(unsigned int i) const {

    if (i >= _fgrids.size())
      return NULL;
    else
      return _fgrids[i];

  }

  std::vector<double> getDomainMin() {
    return _domainMin;
  }

  std::vector<double> getDomainMax() {
    return _domainMax;
  }

  _Tp getCoef(unsigned int i) const {

    if (i >= _fgrids.size())
      return 0.0;
    else
      return _fgrids[i]->getCoef();

  }

  std::vector<FGridContainer<_Tp>*> getFullgrids() {
    return _fgrids;
  }


  int getDim() {

    return _dim;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// VIRTUAL METHODS //////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  /** sets the domain of all the full grids, this is the correct way for extrapolation */
  virtual void setDomainAllFG(GridDomain* gridDomain) const = 0;

  /** create the actual vector for the full grids. <br>
   * This is be different for serial and parallel implementations
   * consider delegating to the current scheme...  */
  virtual void createFullGrids() = 0;

  /** create SGpp grid storage out of the combi grid <br>
   * @return the created grid storage for the SGpp grid
   *
   * */


  ///////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Data Members   //////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
 protected:

  /** pointer to the combi scheme which is set from the constructor argument */
  AbstractCombiScheme<_Tp>* _combischeme;

  /** dimensions of the full grids */
  int _dim;

  /** the coefficients of the full grids*/

  std::vector<FGridContainer<_Tp>*> _fgrids;

  //  std::vector<int> _levels_truncation;

  /** for each dimensions if there are boundary points in that dimension*/
  std::vector<bool> _hasBoundaryPts;

  std::vector<double> _domainMin;

  std::vector<double> _domainMax;

  /**
   * Return the currently used grid stretching type
   */
  Stretching _stretchingType;

  ///////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Protected util methods //////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////

  /** the full grids which appear twice will be deleted.
   * The last instance will be deleted including the coefficient */
  void deleteDuplicate() {

    std::vector<bool> markForDelete(_fgrids.size(), false);
    bool isEqual = false;

    // compare each grid to each grid and check if two are equal
    for (unsigned i = 0; i < _fgrids.size(); i++) {
      for (int j = i + 1; j < _fgrids.size(); j++) {
        FullGrid<double>* fg1 = _fgrids[i]->fg();
        FullGrid<double>* fg2 = _fgrids[j]->fg();
        // test if full grid i is equal
        isEqual = true;

        for (int k = 0; k < _dim; k++) {
          isEqual = (isEqual
                     && (fg1->getLevels()[k] == fg2->getLevels()[k]));
        }

        markForDelete[j] = markForDelete[j] || isEqual;
      }
    }

    // delete the marked full grids
    for (unsigned int i = 0; i < _fgrids.size(); i++) {
      // if this grid was marked then
      if (markForDelete[i]) {
        deleteFullGrid(i);
      }
    }
  }

};

} // namespace combigrid

#endif /* COMBIGRID_HPP_ */
