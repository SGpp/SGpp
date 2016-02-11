// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/combischeme/CombiS_CT.hpp>

#include <vector>

template <typename _Tp>
combigrid::CombiS_CT<_Tp>::CombiS_CT(std::vector<int> levels, int trunc_levels) {
  this->_levels = levels;
  unsigned int dim = (unsigned int)levels.size();
  this->_levels_truncation.resize(dim, trunc_levels);
}

template <typename _Tp>
combigrid::CombiS_CT<_Tp>::CombiS_CT(std::vector<int> levels, std::vector<int> trunc_levels) {
  this->_levels = levels;
  this->_levels_truncation = trunc_levels;
}

template <typename _Tp>
void combigrid::CombiS_CT<_Tp>::initCombiGrid(int in_dim,
                                              std::vector<std::vector<int> >& out_levels_vector,
                                              std::vector<_Tp>& out_coefs) {
  std::vector<int> levels_tmp = this->_levels;
  // the rations for the dimension adaptive case
  std::vector<double> ratio_(in_dim, 1.0);
  // the truncation level, if there is none then l_user = 1 for all dimensions (
  // trapezoidal ) for l_user = 0 (linear)
  std::vector<int> l_user_ = this->_levels_truncation;
  // call the init function
  applyScheme(in_dim, levels_tmp, ratio_, l_user_, out_levels_vector, out_coefs);
}

template <typename _Tp>
void combigrid::CombiS_CT<_Tp>::re_initCombiGrid(int in_dim,
                                                 const std::vector<FGridContainer<_Tp>*> in_grids,
                                                 std::vector<std::vector<int> >& out_levels_vector,
                                                 std::vector<_Tp>& out_coefs) {
  std::vector<std::vector<int> > tmp_levels_vector;
  std::vector<_Tp> tmp_coefs;
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
void combigrid::CombiS_CT<_Tp>::recomputeCoefficients(
    int in_dim, std::vector<FGridContainer<_Tp>*>& out_fgrids) {
  { std::cout << " combiS_CT scheme -> recomputeCoefficients has been invoked \n"; }
}  // do nothing?

/** utility functions for the current implementation of the combination scheme
 **/
template <typename _Tp>
void combigrid::CombiS_CT<_Tp>::getTrapezoidsums(
    std::vector<int>& v, size_t dim, int sum, std::vector<double>& ratio_,
    std::vector<int>& l_user_, std::vector<std::vector<int> >& out_levels_vector) {
  /* Takes recursively every possible combination of numbers which add up to sum
   * creating a linear boundary grid for each one
   * The levels of the full grids must be greater than l_user*/
  // code from Aliz
  if (dim == 1) {
    int tmp = static_cast<int>(sum / ratio_[v.size()]) + l_user_[v.size()];
    v.push_back(tmp);
    // add v to the level vectors
    out_levels_vector.push_back(v);
    v.pop_back();
  } else {
    for (int i = 0; i <= sum; i++) {
      int tmp = static_cast<int>(i / ratio_[v.size()]) + l_user_[v.size()];
      v.push_back(tmp);
      getTrapezoidsums(v, dim - 1, sum - i, ratio_, l_user_, out_levels_vector);
      v.pop_back();
    }
  }
}

template <typename _Tp>
void combigrid::CombiS_CT<_Tp>::applyScheme(int in_dim, const std::vector<int>& in_levels,
                                            std::vector<double>& in_ratio_,
                                            std::vector<int>& in_l_user_,
                                            std::vector<std::vector<int> >& out_levels_vector,
                                            std::vector<_Tp>& out_coefs) {
  std::vector<int> v(0);

  // the ratio for dimension adaptivity
  in_ratio_.resize(in_dim, 1.0);

  // get the maximum level and get the
  int n = in_levels[0], max = in_l_user_[0];

  for (int i = 1; i < in_dim; i++) {
    if (in_levels[i] > n) n = in_levels[i];

    if (in_l_user_[i] > max) max = in_l_user_[i];
  }

  for (int i = 0; i < in_dim; i++) {
    in_ratio_[i] = static_cast<double>(n) / static_cast<double>(in_levels[i]);
  }

  n = n - max;
  // add the different full grids

  _Tp combi = 0.0;

  for (int d = 0; d < in_dim; d++) {
    combi = (_Tp)combigrid::combination(in_dim - 1, d);
    int oldsize = static_cast<int>(out_levels_vector.size());

    if (d % 2 != 0) combi = -combi;

    // call the recursive function to add the spaces
    this->getTrapezoidsums(v, in_dim, n - d, in_ratio_, in_l_user_, out_levels_vector);

    // add the coefficients
    for (int i = oldsize; i < static_cast<int>(out_levels_vector.size()); i++) {
      out_coefs.push_back(combi);
    }
  }

  // remove the duplicated spaces -- > implemented in CombiSchemeBasis...
  this->removeDuplicates(out_levels_vector, out_coefs);
}

template class combigrid::CombiS_CT<float>;
template class combigrid::CombiS_CT<double>;

// add more declarations at your hearth's will !!!!
