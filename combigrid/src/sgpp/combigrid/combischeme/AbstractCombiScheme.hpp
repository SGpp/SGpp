/*
 * AbstractCombiScheme.hpp
 *
 *  Created on: May 22, 2014
 *      Author: petzko
 */

#ifndef ABSTRACTCOMBISCHEME_HPP_
#define ABSTRACTCOMBISCHEME_HPP_

#include <sgpp/combigrid/domain/CombiGridDomain.hpp>
#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>
#include <sgpp/combigrid/fullgrid/GridContainer.hpp>
#include <stdlib.h>
#include "../utils/combigrid_utils.hpp"

namespace combigrid {

template<typename _Tp>
class AbstractCombiScheme {

 protected:

  std::vector<int> _levels;
  std::vector<int> _levels_truncation;

 public:

  /*
   * Replace the internal levels vector of the scheme with new data
   * @param levels - a vector(list) of integers specifying the max level for each dimension
   *
   */
  void setLevels(std::vector<int> levels) {
    _levels = levels;
  }

  /*setter for the truncation levels vector
   * @param trunc_levels: a vector(list) of intergers specifying the truncation level per dimension
   * */
  void setTruncationLevels(std::vector<int> trunc_levels) {
    _levels_truncation = trunc_levels;
  }

  /**
   * Default destructor
   *
   **/
  virtual ~AbstractCombiScheme() {
    _levels.clear();
    _levels_truncation.clear();
  }

  /** The interface function that all classes inheriting from AbstractCombiScheme should implement.
   * initCombiGrid implements the desired combiGrid scheme, i.e. method of construction of the functional spaces
   * and the corresponding coefficient for each space. Notice that the combischeme DOES NOT own any of the data
   * that is given to it. It simply fills in the "out_levels_vector" and "out_coefs" vector data and returns without
   * keeping an internal copies of it. This is a nice way to enforce decoupling between application logic implementation and
   * data management!!
   *
   *
   * @param in_dim (in) an input parameter specifying the dimension of the problem
   * @param out_levels_vector (out) a vector containing the selected levels by the current combigrid implementation
   * @param out_coefs (out) a vector containing the interpolation coeficients for each level...
   *
   * */

  virtual void initCombiGrid(int in_dim,
                             std::vector<std::vector<int> >& out_levels_vector,
                             std::vector<_Tp>& out_coefs) = 0;

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

  virtual void re_initCombiGrid(int in_dim,
                                const std::vector<FGridContainer<_Tp>*> in_grids,
                                std::vector<std::vector<int> >& out_levels_vector,
                                std::vector<_Tp>& out_coefs) = 0;

  /** Implement this method with desired logic to handle situations when recomputation of the coefficients might be necessary.
   *  Exemplary use cases could be the addition or removal of a fullgrid to/from the combigrid container.
   * @param in_dim - dimension of the problem
   * @param out_fgrids - a vector containing the fullgrids associated with the current combigrid- will be updated with the newly recomputed
   * coefficients
   * */
  virtual void recomputeCoefficients(int in_dim,
                                     std::vector<FGridContainer<_Tp>*>& out_fgrids) = 0;

  /**
   * Evaluate the combi grid at one specified point. Buffered evaluation
   * might be faster than one evaluation point.
   * @param in_fgrids - a list of FGridContainers, to be evaluated onto the in_coords
   * @param in_coords the coordinates to evaluate the combigrid on
   *
   * @return the evaluation result!
   */

  _Tp evalCombiGrid(const std::vector<FGridContainer<_Tp>*>& in_fgrids,
                    const std::vector<double>& in_coords) const {

    _Tp result = 0.0;

    // spawn a parallel region for for loop parallelization - due to the different workload for grid evaluation
    // set the scheduling to dynamic with chunk size 1 in effect implementing first-come first serve
    // workload distribution.
    #pragma omp parallel
    {

      #pragma omp for reduction(+:result)  schedule(guided,1)

      for (unsigned int i = 0; i < in_fgrids.size(); i++) {

        // we evaluate each full grid and multiply with the coefficient, and sum the result up
        if (in_fgrids[i]->isActive()) { //evaluate only if grid is active!
          // copy the coordinates before rescaling them
          std::vector<double> coords_tmp = in_coords;
          _Tp coef = in_fgrids[i]->getCoef();
          result += coef
                    * (_Tp) in_fgrids[i]->fg()->eval(coords_tmp);
        }



      }
    }
    //COMBIGRID_OUT_LEVEL3( 4 , "SerialCombiGrid::eval result=" << result);
    return result;

  }

  /** evaluate the combi grid at one specified point. Buffered evaluation
   * which might be faster than one evaluation point.
   * @param in_fgrids - a list of FGridContainers, to be evaluated onto the in_coords
   * @param in_coords - a list of coordinates to be evaluated
   * @param out_result - vector storing the combigrid evaluation results.
   */
  void evalCombiGrid( const std::vector<FGridContainer<_Tp>*>& in_fgrids,
                      const std::vector<std::vector<double> >& in_coords,
                      std::vector<_Tp>& out_results) const {

    out_results.resize(in_coords.size());

    // just iterate over each point and call the serial evaluation function
    for (int i = 0; i < (int) out_results.size(); i++) {
      out_results[i] = evalCombiGrid( in_fgrids,
                                      in_coords[i]);
    }

  }

  /**Evaluate a grid with the specified index and return the result
   *
   * @param index the index of the FullGrid contained in the combigrid
   * @param in_fgrids - a list of FGridContainers, the index-th element of which is to be evaluated onto the in_coords
   * @param in_coords the coordinate points to evaluate
   */
  _Tp evalSingleGrid(const int index,
                     const std::vector<FGridContainer<_Tp>*>& in_fgrids,
                     const std::vector<double>& in_coords) const {

    std::vector<double> coords_temp = in_coords;
    //do not check if the current grid is active within the combigrid !!!
    return in_fgrids[index]->fg()->eval(coords_temp);

  }

 protected:


  void removeDuplicates(std::vector<std::vector<int> >& out_levels_vector,
                        std::vector<_Tp>& out_coefs) {

    // this is Alize's code 1 to 1
    size_t i = 0;
    size_t j;
    std::vector<std::vector<int> >::iterator gt = out_levels_vector.begin();
    typename std::vector<_Tp>::iterator ct = out_coefs.begin();

    int size_fg = static_cast<int>(out_levels_vector.size());

    while ((int) i < size_fg - 1) {
      j = i + 1;

      while ((int) j < size_fg) {
        if (out_levels_vector[i] == out_levels_vector[j]) {
          out_levels_vector.erase(gt + j);
          out_coefs[i] += out_coefs[j];
          out_coefs.erase(ct + j);
          size_fg--;
        } else
          j++;
      }

      if (out_coefs[i] == 0.0) {
        out_levels_vector.erase(gt + i);
        out_coefs.erase(ct + i);
        size_fg--;

        if (i == 0) {
          gt = out_levels_vector.begin();
          ct = out_coefs.begin();
        }
      } else
        i++;
    }
  }

  std::vector<int> updateScheme(std::vector<std::vector<int> > levelsNew,
                                std::vector<_Tp> coefs_new,
                                std::vector<std::vector<int> >& out_levels_vector,
                                std::vector<_Tp>& out_coefs) {

    removeDuplicates(out_levels_vector, out_coefs);
    std::vector<int> result(0);

    for (unsigned int i = 0; i < levelsNew.size(); i++) {
      bool found = false;
      int found_ind = -1;

      for (unsigned int j = 0; j < out_levels_vector.size(); j++) {
        bool equal = true;

        for (unsigned int k = 0; k < out_levels_vector[j].size(); ++k) {
          if (levelsNew[i][k] != out_levels_vector[j][k]) {
            equal = false;
            break;
          }

        }

        if (equal == true) {
          found = true;
          found_ind = j;
          result.push_back(j);
          break;
        }
      }

      if (found) {
        out_coefs[found_ind] += coefs_new[i];
      } else {
        out_levels_vector.push_back((std::vector<int>) levelsNew[i]);
        out_coefs.push_back(coefs_new[i]);
        result.push_back(static_cast<int>(out_coefs.size()) - 1);
      }
    }

    return result;
  }

};
}
#endif /* ABSTRACTCOMBISCHEME_HPP_ */
