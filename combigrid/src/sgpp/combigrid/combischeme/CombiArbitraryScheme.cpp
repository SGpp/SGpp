// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp>

#include <vector>

// using namespace combigrid;
template <typename _Tp>
combigrid::CombiArbitraryScheme<_Tp>::CombiArbitraryScheme(
    std::vector<std::vector<int> > levels_vector) {
  this->_levels = levels_vector[0];
  _levels_vector = levels_vector;
}

template <typename _Tp>
void combigrid::CombiArbitraryScheme<_Tp>::initCombiGrid(
    int in_dim, std::vector<std::vector<int> >& out_levels_vector, std::vector<_Tp>& out_coefs) {
  // acta simplitica!
  CombigridLevelVector buffer = CombigridLevelVector::getCombiLevels(_levels_vector);
  out_levels_vector = buffer.getLevelVec();
  std::vector<double> tmp_coef = buffer.getCoef();
  out_coefs.clear();

  for (unsigned int i = 0; i < tmp_coef.size(); i++) out_coefs.push_back((_Tp)tmp_coef[i]);
}
template <typename _Tp>
void combigrid::CombiArbitraryScheme<_Tp>::re_initCombiGrid(
    int in_dim, const std::vector<FGridContainer<_Tp>*> in_grids,
    std::vector<std::vector<int> >& out_levels_vector, std::vector<_Tp>& out_coefs) {
  std::vector<_Tp> tmp_coefs;
  std::vector<std::vector<int> > tmp_levels_vector;
  initCombiGrid(in_dim, tmp_levels_vector, tmp_coefs);

  // examine what vectors do we have in the in_grids vector....
  // 1) first deactivate all existing grids

  for (unsigned int i = 0; i < in_grids.size(); i++) in_grids[i]->deactivateGrid();

  // 2) if the scheme attempts to create an already existing grid, change the
  // coeffs and activate it.
  for (unsigned int i = 0; i < in_grids.size(); i++) {
    unsigned int j = 0;

    while (j < tmp_levels_vector.size()) {
      // if the grid vector i's levels vector == temp_levels_vector[j] and the
      // coeffs are the same
      if (in_grids[i]->getFGLevels() == tmp_levels_vector[j]) {
        // if true leave the grid as activated..change its coefficient
        // and remove the record from the list of grids to be created.
        in_grids[i]->setCoef(tmp_coefs[j]);
        in_grids[i]->activateGrid();
        unsigned int nr = (unsigned int)tmp_levels_vector.size();
        tmp_levels_vector[j] = tmp_levels_vector[nr - 1];
        tmp_levels_vector.resize(nr - 1);
        tmp_coefs[j] = tmp_coefs[nr - 1];
        tmp_coefs.resize(nr - 1);
      } else {
        j++;
      }
    }
  }

  out_levels_vector = tmp_levels_vector;
  out_coefs = tmp_coefs;
  // now initialize the combigrid as if the in_grids vector were empty..
}

template <typename _Tp>
void combigrid::CombiArbitraryScheme<_Tp>::recomputeCoefficients(
    int in_dim, std::vector<FGridContainer<_Tp>*>& out_fgrids) {
  std::cout << "TODO: IMPLEMENT CombiArbitraryScheme's recomputeCoefficients "
               "method \n";
}

template class combigrid::CombiArbitraryScheme<float>;
template class combigrid::CombiArbitraryScheme<double>;
