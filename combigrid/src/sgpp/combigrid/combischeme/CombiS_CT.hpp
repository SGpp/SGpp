/*
 * CombiS_CT_demo.hpp
 *
 *  Created on: May 22, 2014
 *      Author: petzko
 */

#ifndef COMBIS_CT_HPP_
#define COMBIS_CT_HPP_

#include <sgpp/combigrid/combischeme/AbstractCombiScheme.hpp>

namespace combigrid {

template<typename _Tp>
class CombiS_CT: public AbstractCombiScheme<_Tp> {

 public:

  CombiS_CT(std::vector<int> levels, int trunc_levels = 1);

  CombiS_CT(std::vector<int> levels, std::vector<int>  trunc_levels );

  void initCombiGrid(int in_dim,
                     std::vector<std::vector<int> >& out_levels_vector,
                     std::vector<_Tp>& out_coefs);

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

  ~CombiS_CT() {
  }

 private:
  /** utility functions for the current implementation of the combination scheme
   **/
  void getTrapezoidsums(std::vector<int>& v, size_t dim, int sum,
                        std::vector<double>& ratio_, std::vector<int>& l_user_,
                        std::vector<std::vector<int> >& out_levels_vector);

  void applyScheme(int in_dim, const std::vector<int>& in_levels,
                   std::vector<double>& in_ratio_, std::vector<int>& in_l_user_,
                   std::vector<std::vector<int> >& out_levels_vector,
                   std::vector<_Tp>& out_coefs);

};

}

#endif /* COMBIS_CT_HPP_ */
