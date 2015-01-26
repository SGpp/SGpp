/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "base/grid/generation/functors/ANOVACoarseningFunctor.hpp"
#include <algorithm>
#include <vector>

namespace sg {
  namespace base {

    ANOVACoarseningFunctor::ANOVACoarseningFunctor(DataVector* alpha, size_t
        /*removements_num*/, double threshold, GridStorage* storage) : alpha(alpha),
      threshold(threshold) {
      removements_num = 0;
      int num_anova_components = 1 << storage->dim();
      anova_variances = std::vector<tANOVAValues>(num_anova_components);
      anova_variances_pointers = new tANOVAValues*[num_anova_components];
      std::vector<tANOVAValues>::iterator it;
      int i = 0;

      for (it = anova_variances.begin(); it != anova_variances.end(); it++) {
        it->component_index = i++;
        it->value = 0.0;
        it->points_num = 0;
      }

      GridStorage::index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();
      double L2_value;
      double total_L2_value = 0.0;
      size_t seq;

      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
        index = *(iter->first);
        seq = (*storage)[&index];
        int level_sum = index.getLevelSum();

        if (level_sum > static_cast<int>(storage->dim())) {
          //Ignore the ANOVA Component f0
          double dim = static_cast<double>(storage->dim());
          L2_value = fabs(alpha->get(seq))
                     * pow(2.0, static_cast<double>(dim - level_sum) / 2)
                     / pow(3.0, dim / 2.0);
          int anove_component = getANOVAComponentIndex(index);
          anova_variances[anove_component].value += L2_value;
          anova_variances[anove_component].points_num += 1;
          total_L2_value += L2_value;
        }
      }

      std::sort(anova_variances.begin(), anova_variances.end(), sorter);
      double cut_off_variance = threshold * total_L2_value;

      double cumul_variance = 0.0;
      i = 0;

      for (it = anova_variances.begin(); it != anova_variances.end(); it++) {
        if (cumul_variance <= cut_off_variance) {
          it->refinement_value = 1.0;
        } else {
          it->refinement_value = 0.0;
          removements_num += it->points_num;
        }

        anova_variances_pointers[i++] = &(*it);
        cumul_variance += it->value;
      }


    }

    int ANOVACoarseningFunctor::getANOVAComponentIndex(GridStorage::index_type& index) {
      int comp_index = 0;

      for (size_t d = 0; d < index.dim(); d++) {
        int level = index.getLevel(d);
        int a = level > 1 ? 1 : 0;
        comp_index += (a) * (1 << d);
      }

      return comp_index;
    }

    ANOVACoarseningFunctor::~ANOVACoarseningFunctor() {
      //delete &anova_variances;
      delete [] anova_variances_pointers;
    }


    double ANOVACoarseningFunctor::operator()(GridStorage* storage, size_t seq) {

      return anova_variances_pointers[getANOVAComponentIndex(*(*storage)[seq])]->refinement_value;
    }

    double ANOVACoarseningFunctor::start() {
      return 0.5;
    }

    size_t ANOVACoarseningFunctor::getRemovementsNum() {
      return this->removements_num;
    }

    double ANOVACoarseningFunctor::getCoarseningThreshold() {
      return 0.5;
    }


  } /* namespace base */
} /* namespace sg */
