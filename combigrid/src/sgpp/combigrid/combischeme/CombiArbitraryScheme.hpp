// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIARBITRARYSCHEME_HPP_
#define COMBIARBITRARYSCHEME_HPP_

#include "CombiSchemeBasis.hpp"
#include <sgpp/combigrid/utils/CombigridLevelVector.hpp>

namespace combigrid {
  /** Combischeme created with from an arbitrary active set of full grids. Missing
   * subgrids for the creation of a valid combi solution are automatically added.
   */
  class CombiArbitraryScheme: public combigrid::CombiSchemeBasis {
    public:
      /** Initializes the arbitrary scheme from a set of active set level vectors*/
      CombiArbitraryScheme(std::vector<std::vector<int> > level_vectors);
      /** Initializes the arbitrary scheme from a CombigridLevelVector as acitve set.*/
      CombiArbitraryScheme(CombigridLevelVector input);

  };

}

#endif /* COMBIARBITRARYSCHEME_HPP_ */