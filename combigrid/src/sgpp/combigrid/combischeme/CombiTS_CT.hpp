/*
 * CombiTS_CT.hpp
 *
 *  Created on: 6 Jun 2014
 *      Author: kenny
 */

#ifndef COMBITS_CT_HPP_
#define COMBITS_CT_HPP_

#include <sgpp/combigrid/combischeme/AbstractCombiScheme.hpp>

namespace combigrid {

template<typename _Tp>
class CombiTS_CT: public AbstractCombiScheme<_Tp> {

 private:

  // specifies the minimal hierarchical level for each dimension. the maximal hierarchical level information is
  // stored in the parents' levels vector
  std::vector<int> _levels_small;
  // a flags vector. specifies if the combischeme should be considered for each particular dimension
  std::vector<bool> _makeCombiInDimension;

 public:

  /** Ctor for the TS scheme where the user specifies the max and the min allowed levels
   * for each dimension
   *
   * @param minlevels the min levels
   * @param maxlevels the max levels
   */
  CombiTS_CT(std::vector<int> minlevels, std::vector<int> maxlevels);

  /** Ctor for cases when in specific dimensions no combi should be done
   * @param in_levels the level vector for the dimension adaptive case
   * @param makeCombiInDimension
   * */
  CombiTS_CT(const std::vector<int>& in_levels,
             const std::vector<bool>& makeCombiInDimension);

  /** Simple constructor taking only the levels vector.
   * @param in_max_levels the level vector for the dimension adaptive case
   *
   */
  CombiTS_CT(std::vector<int> in_max_levels);

  /**
   * Set which dimensions should the combischeme be applied to. Default is all dimensions
   * @param makeInDimension <= a flag vector of the size the dimensionality of the problem,
   * with true/false value for each dimension specifying whether or not combischeme should be used for this dim.
   *
   */
  void setActiveDimension(std::vector<bool> makeInDimension) {
    this->_makeCombiInDimension = makeInDimension;
  }

  /** The interface function that all classes inheriting from AbstractCombiScheme should implement.
   * initCombiGrid implements the desired combiGrid scheme, i.e. method of construction of the functional spaces
   * and the corresponding coefficient for each space. Notice that the combischeme DOES NOT own any of the data
   * that is given to it. It simply fills in the "out_levels_vector" and "out_coefs" vector data and returns without
   * keeping an internal copies of it. This is a nice way to enforce decoupling between application logic implementation and
   * data management!!
   *
   *
   * @param in_dim (in) an input parameter specifying the dimensionality of the problem
   * @param out_levels_vector (out) a vector containing the selected levels by the current combigrid implementation
   * @param out_coefs (out) a vector containing the interpolation coefficients for each level...
   *
   * */

  void initCombiGrid(int in_dim,
                     std::vector<std::vector<int> >& out_levels_vector,
                     std::vector<_Tp>& out_coefs);

  /**
   *  re_initCombigrid - a method similar to the initCombiGrid, but this one allows data reusability.
   *  The combischeme will examine the already existing full grids and as a result decide which of the already
   *  existing grids to activate, which to deactivate and which to create anew! The output out_levels_vector and out_coefs vectors
   *  contain the information ONLY for the ADDITIONAL grids that need to be created.
   *
   * @param in_dim (in) an input parameter specifying the dimension of the problem
   * @param in_grids the already existing fullgrids. The scheme uses the information contained inside "in_grids" to decide which grids to activate and which to deactivate
   * @param out_levels_vector (out) a vector containing the selected levels by the current combigrid implementation
   * @param out_coefs (out) a vector containing the interpolation coefficients for each level...
   *
   * */
  void re_initCombiGrid(int in_dim,
                        const std::vector<FGridContainer<_Tp>*> in_grids,
                        std::vector<std::vector<int> >& out_levels_vector,
                        std::vector<_Tp>& out_coefs);

  /** Implement this method with desired logic to handle situations when recomputation of the coefficients might be necessary.
   *  Exemplary use cases could be the addition or removal of a fullgrid to/from the combigrid container.
   * @param in_dim - dimension of the problem
   * @param out_fgrids - a vector containing the fullgrids associated with the current combigrid- will be updated with the newly recomputed
   * coefficients
   * */
  void recomputeCoefficients(int in_dim,
                             std::vector<FGridContainer<_Tp>*>& out_fgrids);

  virtual ~CombiTS_CT() {

    // not strictly necessary since memory will automatically be freed once the class is out of scope... but
    // one can never be too sure :)
    _levels_small.clear();
    this->_levels.clear();
    _makeCombiInDimension.clear();

  }

};

}

#endif /* COMBITS_CT_HPP_ */
